"""Display copy of the raw voltage for the raster viewer's underlay: already filtered and
referenced, in small tiles, at two time resolutions, so that a view downloads little and the
browser has nothing to compute.

    python make_display.py meta    --work D:/temp/fd016_disp [--raw Y:/temp/fd016_raw_v1]
    python make_display.py convert --work D:/temp/fd016_disp [--keys KEYS.json] [--local]
                                   [--t0 1120 --t1 1130] [--probes 0] [--workers 13]
    python make_display.py refs    --work D:/temp/fd016_disp
    python make_display.py sync    --work D:/temp/fd016_disp --keys KEYS.json

It goes with the copy of the recorded samples written by make_raw.py (prefix v1/ in the
bucket), whose manifest supplies the probes, dead sites, scaling, clock and motion files; this
copy lives under PREFIX next to it and adds a manifest of its own.

What is stored. The viewer's display processing, done here once, from the original recording
(not from the 2-bits-dropped copy):
  1. per site, one second at a time with MARGIN samples (10 ms) on each side: subtract the mean
     of that stretch, high-pass (Butterworth, order HP_ORDER at HP_HZ) forward and backward
  2. dead sites are set to 0; per sample, the median over the shank's other sites is subtracted
  3. microvolts, divided by STEP_UV, rounded, clipped to +-127: one signed byte per value
Level 0 keeps every sample, in tiles of 50 ms. Level 1 keeps, for every block of 8 samples,
the sample of largest magnitude (the viewer draws the sample of largest magnitude in each
pixel, so this is what a pixel that covers whole blocks would show), in tiles of 0.5 s.
Per shank and level one Zarr (version 2) array, 384 sites x values, int8, zlib; no filter, so
zarr-python reads it as it is and a browser decodes a chunk with DecompressionStream("deflate").
Tile T of a level holds values T*tile to T*tile + tile - 1, counted from sample 0 of the
probe's recording (value j of level 1 is samples 8j to 8j + 7).

Nothing is staged on disk unless --local is given (used for a short stretch to test the viewer
with serve.py --disp): tiles go straight to the bucket. progress/p<probe>_sh<shank>.jsonl in the
work folder records every second done, with each tile's size, MD5 and sum. `sync` lists the
bucket, compares name, size and MD5, converts again whatever is missing or different, and fails
while anything of the manifest's range is absent. Only when everything is there does it upload
the small files, the manifest last: the viewer takes the manifest as the sign that the copy is
complete.

The compressed bytes depend on the zlib build (the stored values do not), so `sync` decides by
the values: a tile converted again must have the same values as the first time (their MD5,
"vmd5", where the record has one, otherwise their sum), and if only its bytes differ, the new
size and MD5 are recorded. A tile that is already in the bucket is not replaced by one with
other values: the values of a stored tile must never change under the same PREFIX, because
browsers keep tiles for a year. A second that fails (a read error, say) is reported and left
for the next run; convert and sync exit with an error when anything is missing. Do not run
convert and sync at the same time.
"""
import argparse
import hashlib
import json
import sys
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from cbin import Cbin
from make_raw import CACHE_CHUNK, CACHE_SMALL, s3_client
from make_raw_proto import FS, HP_HZ, HP_ORDER, MARGIN, highpass_sos, raw_path

PREFIX = "d1/"                  # folder in the bucket; a changed format gets a new one
STEP_UV = 6.0                   # microvolts per stored count
CLIP = 127
ZLIB_LEVEL = 6
LEVELS = [dict(level=0, block=1, tile=1500),      # values per tile; a tile is tile * block samples
          dict(level=1, block=8, tile=1875)]
NCH = 384
TASK_S = 10                     # seconds of one shank handled by one task (they share the reads)


def tile_path(p, sh, level, T):
    return f"p{p}/sh{sh}/L{level}/0.{T}"


def process_second(get, c, n_chunks, dead, uv_per_count, sos=None):
    """Second c (samples c*FS .. c*FS + FS - 1) of one shank -> float64 microvolts (384, FS).

    get(k) returns chunk k of the recording as int16 (FS, >= 384 channels).
    """
    from scipy import signal
    parts, lo = [], 0
    if c > 0:
        parts.append(get(c - 1)[-MARGIN:, :NCH]); lo = MARGIN
    parts.append(get(c)[:, :NCH])
    n = parts[-1].shape[0]                       # the recording's last second may be short
    if c + 1 < n_chunks:
        parts.append(get(c + 1)[:MARGIN, :NCH])
    x = np.concatenate(parts).T.astype(np.float64)
    x -= x.mean(axis=1, keepdims=True)
    sos = highpass_sos() if sos is None else sos
    y = signal.sosfilt(sos, x, axis=1)
    y = signal.sosfilt(sos, y[:, ::-1], axis=1)[:, ::-1]
    good = np.setdiff1d(np.arange(NCH), dead)
    y[dead] = 0.0
    y[good] -= np.median(y[good], axis=0, keepdims=True)
    out = np.zeros((NCH, FS))
    out[:, :n] = y[:, lo:lo + n] * uv_per_count
    return out


def quantize(y):
    return np.clip(np.rint(y / STEP_UV), -CLIP, CLIP).astype(np.int8)


def block_extreme(y, k):
    """Per block of k samples along the last axis, the sample of largest magnitude (the first
    one if several are equal)."""
    b = y.reshape(y.shape[0], -1, k)
    i = np.abs(b).argmax(axis=2)
    return np.take_along_axis(b, i[..., None], axis=2)[..., 0]


def tiles_of_second(y, c):
    """Float microvolts of second c -> list of (level, T, int8 (384, tile))."""
    out = []
    for lv in LEVELS:
        v = quantize(y if lv["block"] == 1 else block_extreme(y, lv["block"]))
        per = v.shape[1] // lv["tile"]
        assert per * lv["tile"] == v.shape[1]
        for i in range(per):
            out.append((lv["level"], c * per + i, np.ascontiguousarray(v[:, i * lv["tile"]:(i + 1) * lv["tile"]])))
    return out


def encode_tile(v):
    return zlib.compress(v.tobytes(), ZLIB_LEVEL)


def decode_tile(buf, tile):
    a = np.frombuffer(zlib.decompress(buf), np.int8)
    assert a.size == NCH * tile
    return a.reshape(NCH, tile)


# ---------------------------------------------------------------------------- meta
def second_span(pr, raw_tile):
    """First and last second (1 s chunk of the recording) that hold samples of the v1 copy's tiles."""
    return (pr["tile_first"] * raw_tile) // FS, ((pr["tile_last"] + 1) * raw_tile - 1) // FS


def do_meta(a):
    import shutil
    work = Path(a.work)
    (work / "v1").mkdir(parents=True, exist_ok=True)
    raw = json.load(open(Path(a.raw) / "manifest.json"))
    for f in ["manifest.json"] + [pr["clock"] for pr in raw["probes"]]:      # kept with the progress files
        shutil.copyfile(Path(a.raw) / f, work / "v1" / f)
    assert raw["version"] == 2 and raw["fs"] == FS and raw["n_channels"] == NCH and raw["margin"] == MARGIN
    assert np.allclose(raw["highpass"]["sos"], highpass_sos())
    for lv in LEVELS:
        assert FS % (lv["tile"] * lv["block"]) == 0, "tiles must divide a second"
    probes = []
    for pr in raw["probes"]:
        cb = Cbin(raw_path(a.session, pr["probe"], 0))
        assert np.all(np.diff(cb.bounds[:-1]) == FS), "the .cbin chunks are expected to be 1 s"
        assert cb.n_samples == pr["n_samples"]
        lvs = []
        c0, c1 = second_span(pr, raw["tile"])
        assert c1 < cb.n_chunks
        for lv in LEVELS:
            per = FS // (lv["tile"] * lv["block"])           # whole seconds are stored
            lvs.append(dict(tile_first=int(c0 * per), tile_last=int((c1 + 1) * per - 1)))
            for s in pr["shanks"]:
                d = work / f"p{pr['probe']}/sh{s['shank']}/L{lv['level']}"
                d.mkdir(parents=True, exist_ok=True)
                n_values = -(-pr["n_samples"] // (lv["tile"] * lv["block"])) * lv["tile"]
                json.dump(dict(zarr_format=2, shape=[NCH, n_values], chunks=[NCH, lv["tile"]], dtype="|i1",
                               compressor=dict(id="zlib", level=ZLIB_LEVEL), fill_value=0, order="C",
                               filters=None, dimension_separator="."), open(d / ".zarray", "w"), indent=1)
                json.dump(dict(step_uv=STEP_UV, block=lv["block"], fs=FS,
                               note="microvolts = value * step_uv; value j is the sample of largest "
                                    "magnitude among samples j*block .. j*block + block - 1"),
                          open(d / ".zattrs", "w"), indent=1)
        probes.append(dict(probe=pr["probe"], name=pr["name"], n_samples=pr["n_samples"],
                           second_first=int(c0), second_last=int(c1), levels=lvs,
                           shanks=[dict(shank=s["shank"]) for s in pr["shanks"]]))
    man = dict(version=1, kind="display",
               description="Raw voltage of FD_016 2026-11-13/1 as the raster viewer displays it "
                           "(high-passed, median-referenced, in steps of STEP_UV); see "
                           "raw_underlay/make_display.py in github.com/SteinmetzLab/rasterViewer",
               fs=FS, n_channels=NCH, step_uv=STEP_UV, clip=CLIP, margin=MARGIN,
               highpass=dict(hz=HP_HZ, order=HP_ORDER),
               levels=[dict(lv, path="p{probe}/sh{shank}/L%d/0.{tile}" % lv["level"]) for lv in LEVELS],
               probes=probes, test=[])
    json.dump(man, open(work / "manifest.json", "w"), indent=1)
    print(f"manifest for {len(probes)} probes written to {work}")


def load(work):
    work = Path(work)
    man = json.load(open(work / "manifest.json"))
    raw_dir = work / "v1"
    raw = json.load(open(raw_dir / "manifest.json"))
    return man, raw, raw_dir


# ---------------------------------------------------------------------------- convert
_W = {}


def check_zlib(man):
    """The bytes of a tile depend on the zlib build; say so when it is not the one of `refs`."""
    was = man.get("zlib_runtime")
    if was and was != zlib.ZLIB_RUNTIME_VERSION:
        print(f"note: this Python's zlib is {zlib.ZLIB_RUNTIME_VERSION}, the manifest was written with {was}; "
              f"tiles made now have other bytes (the same values)", flush=True)


def _init(session, work, keys, local):
    _W.update(session=session, work=Path(work), local=local, cb={}, s3=None, pool=None, raw=load(work)[1])
    if keys:
        _W["s3"], _W["bucket"] = s3_client(keys)
        _W["pool"] = ThreadPoolExecutor(10)


def _put(key, body):
    try:
        _W["s3"].put_object(Bucket=_W["bucket"], Key=key, Body=body,
                            ContentType="application/octet-stream", CacheControl=CACHE_CHUNK)
        return True
    except Exception as e:                                   # `sync` converts it again
        sys.stderr.write(f"upload of {key} failed: {type(e).__name__}: {e}\n")
        return False


def _task(task):
    """Some consecutive seconds of one shank -> their tiles, uploaded (and written with --local).

    only: None, or {(level, T): expected} for the tiles to store, the others being computed
    but left alone. expected is None, or the (vmd5 or None, sum) the tile must have; a tile
    with other values is not stored and comes back with mismatch=True.
    """
    p, sh, seconds, only = task
    cb = _W["cb"].get((p, sh))
    if cb is None:
        cb = _W["cb"][(p, sh)] = Cbin(raw_path(_W["session"], p, sh))
    s = _W["raw"]["probes"][p]["shanks"][sh]
    cache = {}

    def get(k):
        if k not in cache:
            cache[k] = cb.chunk(k)
            for old in [q for q in cache if q < k - 2]:
                del cache[old]
        return cache[k]

    res = []
    for c in seconds:
        try:
            y = process_second(get, c, cb.n_chunks, s["dead"], s["uv_per_count"])
        except Exception as e:                               # reported, and left for the next run
            res.append(dict(c=int(c), error=f"{type(e).__name__}: {e}"))
            continue
        recs, puts = [], []
        for level, T, v in tiles_of_second(y, c):
            buf = encode_tile(v)
            r = dict(L=level, T=int(T), n=len(buf), md5=hashlib.md5(buf).hexdigest(),
                     sum=int(v.sum(dtype=np.int64)), vmd5=hashlib.md5(v.tobytes()).hexdigest(), up=False)
            recs.append(r)
            if only is not None:
                if (level, int(T)) not in only:
                    puts.append(None)
                    continue
                exp = only[(level, int(T))]
                if exp is not None and (r["sum"] != exp[1] or (exp[0] and r["vmd5"] != exp[0])):
                    r["mismatch"] = True
                    puts.append(None)
                    continue
            rel = tile_path(p, sh, level, T)
            if _W["local"]:
                (_W["work"] / rel).write_bytes(buf)
            puts.append(_W["pool"].submit(_put, PREFIX + rel, buf) if _W["s3"] else None)
        for r, f in zip(recs, puts):
            r["up"] = bool(f and f.result())
        res.append(dict(c=int(c), tiles=recs))
    return p, sh, res


def read_progress(work, p, sh):
    """-> {second: record}; a later record of the same second replaces an earlier one."""
    f = Path(work) / "progress" / f"p{p}_sh{sh}.jsonl"
    done = {}
    if f.exists():
        for line in open(f):
            line = line.strip()
            if line:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:                 # a line cut short by a crash
                    continue
                done[r["c"]] = r
    return done


def second_range(man, raw, raw_dir, pr, t0=None, t1=None):
    """First and last second (chunk of the recording) to convert for a probe."""
    lo, hi = pr["second_first"], pr["second_last"]
    if t0 is not None:
        k = np.fromfile(raw_dir / raw["probes"][pr["probe"]]["clock"], "<f8").reshape(-1, 2)
        lo = max(lo, int(np.interp(t0, k[:, 1], k[:, 0]) // FS))
        hi = min(hi, int(np.interp(t1, k[:, 1], k[:, 0]) // FS))
    return lo, hi


def append_progress(work, p, sh, results):
    pf = Path(work) / "progress" / f"p{p}_sh{sh}.jsonl"
    cut = pf.exists() and pf.stat().st_size > 0 and not pf.read_bytes().endswith(b"\n")
    with open(pf, "a") as log:
        if cut:                                              # a line cut short by a crash
            log.write("\n")
        for r in results:
            log.write(json.dumps(r) + "\n")


def do_convert(a):
    import multiprocessing as mp
    work = Path(a.work)
    man, raw, raw_dir = load(work)
    (work / "progress").mkdir(exist_ok=True)
    assert a.keys or a.local, "give --keys (upload), --local (write files), or both"
    t_begin = time.time()
    n_sec = n_bytes = n_fail = n_err = 0
    check_zlib(man)
    with mp.Pool(a.workers, initializer=_init, initargs=(a.session, str(work), a.keys, a.local)) as pool:
        for pr in man["probes"]:
            p = pr["probe"]
            if a.probes is not None and p not in a.probes:
                continue
            lo, hi = second_range(man, raw, raw_dir, pr, a.t0, a.t1)
            for s in pr["shanks"]:
                sh = s["shank"]
                done = read_progress(work, p, sh)
                todo = [c for c in range(lo, hi + 1) if c not in done]
                tasks, run = [], []
                for c in todo:                               # runs of consecutive seconds, TASK_S at most
                    if run and (c != run[-1] + 1 or len(run) == TASK_S):
                        tasks.append((p, sh, run, None)); run = []
                    run.append(c)
                if run:
                    tasks.append((p, sh, run, None))
                if not tasks:
                    print(f"probe {p} shank {sh}: nothing to do", flush=True)
                    continue
                t_sh, b_sh = time.time(), 0
                for _, _, res in pool.imap_unordered(_task, tasks):
                    for r in res:
                        if "error" in r:
                            n_err += 1
                            print(f"probe {p} shank {sh} second {r['c']}: {r['error']}", flush=True)
                    res = [r for r in res if "error" not in r]
                    append_progress(work, p, sh, res)
                    for r in res:
                        n_sec += 1
                        b_sh += sum(t["n"] for t in r["tiles"])
                        n_fail += sum(1 for t in r["tiles"] if a.keys and not t["up"])
                n_bytes += b_sh
                print(f"probe {p} shank {sh}: {len(todo)} s in {time.time() - t_sh:.0f} s, "
                      f"{b_sh / len(todo) / 1e6:.2f} MB per s of recording; total {n_sec} s, "
                      f"{n_bytes / 1e9:.1f} GB, {n_fail} failed uploads, {n_err} seconds failed, "
                      f"{(time.time() - t_begin) / 60:.1f} min", flush=True)
    print(f"convert finished: {n_sec} shank-seconds, {n_bytes / 1e9:.2f} GB, {n_fail} failed uploads, "
          f"{n_err} seconds failed, {(time.time() - t_begin) / 60:.1f} min", flush=True)
    return n_fail + n_err                      # nonzero: run sync (failed uploads) or convert again


# ---------------------------------------------------------------------------- refs
def do_refs(a):
    """One check sum per shank and level (for the viewer's ?rawtest=1) and the sizes."""
    work = Path(a.work)
    man, raw, raw_dir = load(work)
    n_tiles = [0] * len(LEVELS)
    n_bytes = [0] * len(LEVELS)
    for pr in man["probes"]:
        p = pr["probe"]
        k = np.fromfile(raw_dir / raw["probes"][p]["clock"], "<f8").reshape(-1, 2)
        smp = int(np.interp(a.t_ref, k[:, 1], k[:, 0]))
        for s in pr["shanks"]:
            done = read_progress(work, p, s["shank"])
            recs = {(t["L"], t["T"]): t for r in done.values() for t in r["tiles"]}
            for t in recs.values():
                n_tiles[t["L"]] += 1
                n_bytes[t["L"]] += t["n"]
            s["check"] = []
            for lv in LEVELS:
                T = smp // (lv["tile"] * lv["block"])
                if (lv["level"], T) in recs:
                    s["check"].append(dict(level=lv["level"], tile=T, sum=recs[(lv["level"], T)]["sum"]))
    for lv, nt, nb in zip(man["levels"], n_tiles, n_bytes):
        if nt:
            lv.update(tiles_written=nt, stored_bytes=nb, bytes_per_s=nb / (nt * lv["tile"] * lv["block"] / FS))
    man.setdefault("zlib_runtime", zlib.ZLIB_RUNTIME_VERSION)
    json.dump(man, open(work / "manifest.json", "w"), indent=1)
    print("; ".join(f"level {lv['level']}: {nt} tiles, {nb / 1e9:.2f} GB" for lv, nt, nb in zip(LEVELS, n_tiles, n_bytes)))


# ---------------------------------------------------------------------------- sync
def do_sync(a):
    import multiprocessing as mp
    work = Path(a.work)
    man, raw, raw_dir = load(work)
    s3, bucket = s3_client(a.keys)
    have = {}
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=PREFIX):
        for o in page.get("Contents", []):
            have[o["Key"]] = (o["Size"], o["ETag"].strip('"'))
    print(f"{len(have)} objects in the bucket under {PREFIX}", flush=True)
    check_zlib(man)
    expected, absent, redo, records = set(), 0, {}, {}
    for pr in man["probes"]:
        p = pr["probe"]
        lo, hi = second_range(man, raw, raw_dir, pr)
        for s in pr["shanks"]:
            sh = s["shank"]
            done = read_progress(work, p, sh)
            miss = [c for c in range(lo, hi + 1) if c not in done]
            if miss:
                absent += len(miss)
                print(f"  probe {p} shank {sh}: {len(miss)} s not converted yet (first {miss[0]}, last {miss[-1]})")
            for c, r in done.items():
                marked = False
                for t in r["tiles"]:
                    key = PREFIX + tile_path(p, sh, t["L"], t["T"])
                    expected.add(key)
                    if have.get(key) != (t["n"], t["md5"]):
                        # in the bucket with other bytes: its values must stay; absent: any will do
                        redo.setdefault((p, sh, c), {})[(t["L"], t["T"])] = \
                            (t.get("vmd5"), t["sum"]) if key in have else None
                        records[(p, sh, c)] = r
                    elif not t["up"]:                        # there, though the upload had reported a failure
                        t["up"] = marked = True
                if marked and (p, sh, c) not in redo:
                    append_progress(work, p, sh, [r])
    n_redo = sum(len(v) for v in redo.values())
    print(f"{len(expected)} tiles in the progress files, {n_redo} to convert again and upload; "
          f"{absent} shank-seconds of the manifest's range are not converted yet", flush=True)
    bad = 0
    if redo:
        tasks = [(p, sh, [c], only) for (p, sh, c), only in redo.items()]
        n_other = 0
        with mp.Pool(min(a.workers, len(tasks)), initializer=_init, initargs=(a.session, str(work), a.keys, False)) as pool:
            for p, sh, res in pool.imap_unordered(_task, tasks):
                for r in res:
                    if "error" in r:
                        bad += len(redo[(p, sh, r["c"])])
                        print(f"  probe {p} shank {sh} second {r['c']}: {r['error']}")
                        continue
                    want, new = redo[(p, sh, r["c"])], {(t["L"], t["T"]): t for t in r["tiles"]}
                    rec = records[(p, sh, r["c"])]
                    for t in rec["tiles"]:                   # the second's record, brought up to date
                        k = (t["L"], t["T"])
                        if k not in want:
                            continue
                        if new[k].get("mismatch"):
                            bad += 1
                            print(f"  probe {p} shank {sh} level {k[0]} tile {k[1]}: converted again, the VALUES "
                                  f"differ from those in the bucket; not uploaded, not recorded")
                        elif not new[k]["up"]:
                            bad += 1
                            print(f"  probe {p} shank {sh} level {k[0]} tile {k[1]}: upload failed")
                        else:
                            n_other += new[k]["md5"] != t["md5"]
                            t.update(n=new[k]["n"], md5=new[k]["md5"], sum=new[k]["sum"], vmd5=new[k]["vmd5"], up=True)
                    append_progress(work, p, sh, [rec])
        if n_other:
            print(f"  {n_other} tiles have the same values but other bytes than the first time (another zlib); "
                  f"the new size and MD5 are recorded")
    small = ["manifest.json"] + [f"p{pr['probe']}/sh{s['shank']}/L{lv['level']}/{f}" for pr in man["probes"]
                                 for s in pr["shanks"] for lv in LEVELS for f in (".zarray", ".zattrs")]
    # The manifest tells the viewer that the copy is there, so it goes up last and only when
    # every tile is in the bucket.
    complete = not (bad or absent)
    for rel in small[1:] + small[:1] if complete else []:
        s3.put_object(Bucket=bucket, Key=PREFIX + rel, Body=(work / rel).read_bytes(),
                      ContentType="application/json", CacheControl=CACHE_SMALL)
    expected.update(PREFIX + rel for rel in small)
    stray = sorted(set(have) - expected)
    print(f"converted again and uploaded {n_redo - bad} tiles ({bad} problems); "
          + (f"uploaded {len(small)} small files, the manifest last; " if complete else
             "NOT complete, so the manifest was not uploaded; ")
          + f"{len(stray)} objects in the bucket are not expected{': ' + ', '.join(stray[:5]) if stray else ''}")
    return bad + absent                       # nonzero: the bucket does not hold the whole range yet


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["meta", "convert", "refs", "sync"])
    ap.add_argument("--session", default="Y:/Subjects/FD_016/2026-11-13/1")
    ap.add_argument("--work", required=True, help="folder for the manifest, the progress files and, with --local, the tiles")
    ap.add_argument("--raw", default="Y:/temp/fd016_raw_v1", help="folder of the v1 copy (its manifest and clock files)")
    ap.add_argument("--keys", default=None, help="JSON with the bucket's credentials; convert uploads if given")
    ap.add_argument("--local", action="store_true", help="convert also writes the tiles into the work folder")
    ap.add_argument("--probes", type=int, nargs="*", default=None)
    ap.add_argument("--t0", type=float, default=None, help="convert only this stretch (viewer s)")
    ap.add_argument("--t1", type=float, default=None)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--t-ref", type=float, default=1125.05, help="viewer time of the check tiles")
    a = ap.parse_args()
    assert (a.t0 is None) == (a.t1 is None)
    (Path(a.work) / "progress").mkdir(parents=True, exist_ok=True)
    r = {"meta": do_meta, "convert": do_convert, "refs": do_refs, "sync": do_sync}[a.step](a)
    sys.exit(1 if r else 0)
