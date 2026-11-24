"""Raw voltage of the whole recording for the raster viewer's underlay: all probes, one version.

    python make_raw.py meta    --out Y:/temp/fd016_raw_v1
    python make_raw.py convert --out Y:/temp/fd016_raw_v1 [--keys KEYS.json] [--t0 1095 --t1 1165]
    python make_raw.py refs    --out Y:/temp/fd016_raw_v1
    python make_raw.py sync    --out Y:/temp/fd016_raw_v1 --keys KEYS.json

What is stored. Per shank one Zarr (version 2) array, 384 sites x samples, int16, in chunks of
0.2 s (6000 samples, a "tile"): the recorded samples with their SHIFT lowest bits dropped
(value = sample >> SHIFT; the viewer uses value * 2**SHIFT + (2**SHIFT - 1) / 2). Nothing else
is done to the data: no filtering, no referencing, no motion correction. Each chunk is stored
with a delta filter and zlib, so a browser decodes it with DecompressionStream("deflate") and a
running sum. The chunk files are written directly (not through the zarr package), and
verify_raw.py reads them back through zarr-python as an independent check.

Tiles are numbered from sample 0 of the probe's recording (tile T holds samples T*6000 to
T*6000 + 5999). Only tiles that overlap the synced period (the viewer's time range) are
written, plus one on each side for the display filter; chunks outside are absent, which Zarr
reads as zeros.

`meta` writes what the viewer needs besides the samples:
  manifest.json        probes, shanks, dead sites, scaling, tile ranges, filter coefficients
  clock_p<probe>.f64   pairs (probe sample, time on the viewer's clock), float64; see clock_knots
  motion_p<probe>.f32  the sorter's motion estimate, float32 (shank, time bin, depth window);
                       the viewer shifts the voltage by it, because the sort (and so every
                       unit's depth) is in motion-corrected coordinates
  p<probe>/sh<shank>/.zarray
`convert` writes the chunks (resumable: progress/p<probe>_sh<shank>.jsonl lists the tiles done,
with size, MD5 and the sum of the stored values) and, with --keys, uploads each to the bucket.
`refs` writes reference output of the display processing for a few tiles (test/), which the
viewer's ?rawtest=1 compares with its own JavaScript, and records one chunk sum per shank.
`sync` compares the bucket with the progress files (key, size, MD5) and uploads what is
missing or different, and the small files.

KEYS.json: {"account_id", "access_key_id", "secret_access_key", "bucket"} for Cloudflare R2
(S3 API). It is a secret and is not part of this repository.
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
from make_raw_proto import (FS, HP_HZ, HP_ORDER, MARGIN, PROBES, TILE, clock_pairs,
                            highpass_sos, process_tile, raw_path, sort_dir)

SHIFT = 2                       # lowest bits dropped
PREFIX = "v1/"                  # folder in the bucket; a changed format gets a new one
ZLIB_LEVEL = 4
VIEWER_META = Path(__file__).resolve().parent.parent / "data" / "meta.json"
CACHE_CHUNK = "public, max-age=31536000, immutable"
CACHE_SMALL = "public, max-age=300"


# ---------------------------------------------------------------------------- clock
def clock_knots(session, probe):
    """The map from probe sample to viewer time, as the knots of a piecewise-linear function.

    The sorted spike times were put on the rig's clock by linear interpolation between sync
    edges. The edges' samples are known (<probe>_sync.npy); how the pipeline paired them with
    the rig's edges is not needed, because every sorted spike carries both a sample and a
    time: the time at each edge is extrapolated from two spikes of the interval on its left
    and from two of the interval on its right. Edges with no usable interval next to them are
    dropped, and any spike of such an interval becomes a knot itself. Returns knots (sample,
    time) and a dict of checks: every spike's time is recomputed from the knots.
    """
    ps, pt, skipped = clock_pairs(session, probe)
    edges = np.round(np.load(Path(session) / f"{PROBES[probe]}_sync.npy").ravel() * FS)
    edges = np.unique(edges[(edges > ps[0]) & (edges < ps[-1])])
    ks = np.concatenate([[ps[0]], edges, [ps[-1]]])
    i0 = np.searchsorted(ps, ks[:-1], "left")              # interval j: spikes with
    i1 = np.searchsorted(ps, ks[1:], "right") - 1          # ks[j] <= sample <= ks[j + 1]
    ok = i1 > i0
    ok[ok] = ps[i1[ok]] > ps[i0[ok]]
    slope = np.full(ks.size - 1, np.nan)
    slope[ok] = (pt[i1[ok]] - pt[i0[ok]]) / (ps[i1[ok]] - ps[i0[ok]])
    line = lambda j, x: pt[i0[j]] + slope[j] * (x - ps[i0[j]])
    kt = np.full(ks.size, np.nan)
    kt[0], kt[-1] = pt[0], pt[-1]
    mid = np.arange(1, ks.size - 1)
    vl = np.where(ok[mid - 1], line(mid - 1, ks[mid]), np.nan)
    vr = np.where(ok[mid], line(mid, ks[mid]), np.nan)
    both = np.isfinite(vl) & np.isfinite(vr)
    kt[mid] = np.where(both, 0.5 * (vl + vr), np.where(np.isfinite(vl), vl, vr))
    keep = np.isfinite(kt)
    # an interval without a usable line (no two spikes at different samples) can still hold a
    # spike, which is a point of the map: those become knots too
    extra = np.concatenate([np.arange(i0[j], i1[j] + 1) for j in np.flatnonzero(~ok)] or [[]]).astype(int)
    ks, kt = np.concatenate([ks[keep], ps[extra]]), np.concatenate([kt[keep], pt[extra]])
    ks, first = np.unique(ks, return_index=True)
    kt = kt[first]
    assert np.all(np.diff(ks) > 0) and np.all(np.diff(kt) > 0), "clock knots are not increasing"
    err = (np.interp(ps, ks, kt) - pt) * FS                 # samples
    lin = np.polyfit(ks, kt, 1)
    wob = (kt - np.polyval(lin, ks)) * 1e3
    chk = dict(probe=probe, n_spikes=int(ps.size), units_skipped=int(skipped), n_edges=int(edges.size),
               n_knots=int(ks.size), edges_dropped=int((~keep).sum()), spike_knots=int(np.unique(ps[extra]).size),
               left_right_max_samples=float(np.nanmax(np.abs(vl - vr)) * FS) if both.any() else None,
               spike_err_max_samples=float(np.abs(err).max()),
               spike_err_p999_samples=float(np.quantile(np.abs(err), 0.999)),
               wobble_ms=[float(wob.min()), float(wob.max())],
               sample_range=[float(ks[0]), float(ks[-1])], time_range=[float(kt[0]), float(kt[-1])])
    return ks, kt, chk


# ---------------------------------------------------------------------------- meta
def shank_info(session, probe, sh):
    so = sort_dir(session, probe, sh)
    prov = json.load(open(so / "provenance.json"))
    ko = so / "kilosort4" / "sorter_output"
    pos = np.load(ko / "channel_positions.npy")
    cmap = np.load(ko / "channel_map.npy")
    c = np.arange(384)
    assert np.array_equal(cmap, c), "sorter dropped or reordered channels"
    assert np.array_equal(pos[:, 1], 15.0 * (c // 2)) and np.array_equal(pos[:, 0], pos[c % 2, 0]) \
        and pos[0, 0] < pos[1, 0], "unexpected site geometry"
    assert prov.get("motion_preset") == "dredge" and prov["kilosort_settings"]["nblocks"] == 0
    m = np.load(so / "motion.npz")
    tb, sb, disp = m["temporal_bins_s"], m["spatial_bins_um"], m["displacement"]
    assert disp.shape == (tb.size, sb.size) and np.allclose(np.diff(tb), tb[1] - tb[0]) and np.all(np.isfinite(disp))
    cb = Cbin(raw_path(session, probe, sh))
    assert cb.n_channels == 385 and cb.fs == FS
    assert np.array_equal(cb.bounds[:-1], np.arange(cb.n_chunks) * FS), "chunks are not whole seconds"
    return dict(shank=sh, x_um=[float(pos[0, 0]), float(pos[1, 0])],
                dead=sorted(int(i) for i in prov["destripe"]["dead_channel_indices"]),
                uv_per_count=float(prov["recording"]["uv_per_bit"]), n_samples=cb.n_samples,
                motion=dict(t0_s=float(tb[0]), dt_s=float(tb[1] - tb[0]), depths_um=sb.tolist(), disp=disp))


def zarray(n_samples):
    return {"shape": [384, int(n_samples)], "chunks": [384, TILE], "dtype": "<i2", "fill_value": 0,
            "order": "C", "filters": [{"id": "delta", "dtype": "<i2", "astype": "<i2"}],
            "dimension_separator": ".", "compressor": {"id": "zlib", "level": ZLIB_LEVEL},
            "zarr_format": 2}


def do_meta(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    vm = json.load(open(VIEWER_META))
    ts, te = float(vm["t_start"]), float(vm["t_end"])
    sos = highpass_sos()
    probes, checks = [], []
    for p in range(len(PROBES)):
        ks, kt, chk = clock_knots(a.session, p)
        checks.append(chk)
        print(json.dumps(chk))
        assert chk["spike_err_max_samples"] < 0.05, "the knots do not reproduce the spike times"
        np.stack([ks, kt], axis=1).astype("<f8").tofile(out / f"clock_p{p}.f64")
        shanks = [shank_info(a.session, p, sh) for sh in range(4)]
        n = shanks[0]["n_samples"]
        mo = shanks[0]["motion"]
        for s in shanks:
            assert s["n_samples"] == n and s["motion"]["disp"].shape == mo["disp"].shape
            assert s["motion"]["t0_s"] == mo["t0_s"] and s["motion"]["dt_s"] == mo["dt_s"]
            assert s["motion"]["depths_um"] == mo["depths_um"] and s["x_um"] == shanks[0]["x_um"]
        np.stack([s["motion"]["disp"] for s in shanks]).astype("<f4").tofile(out / f"motion_p{p}.f32")
        # samples at the ends of the viewer's time range (the map extended linearly)
        s_a = ks[0] + (ts - kt[0]) * (ks[1] - ks[0]) / (kt[1] - kt[0])
        s_b = ks[-1] + (te - kt[-1]) * (ks[-1] - ks[-2]) / (kt[-1] - kt[-2])
        n_tiles = -(-n // TILE)
        t_first = max(0, int(s_a // TILE) - 1)
        t_last = min(n_tiles - 1, int(s_b // TILE) + 1)
        for s in shanks:
            d = out / f"p{p}" / f"sh{s['shank']}"
            d.mkdir(parents=True, exist_ok=True)
            json.dump(zarray(n), open(d / ".zarray", "w"), indent=1)
        probes.append(dict(
            probe=p, name=PROBES[p], n_samples=int(n), clock=f"clock_p{p}.f64", n_knots=int(ks.size),
            tile_first=t_first, tile_last=t_last, x_um=shanks[0]["x_um"],
            motion=dict(file=f"motion_p{p}.f32", shape=[4, *mo["disp"].shape], t0_s=mo["t0_s"],
                        dt_s=mo["dt_s"], depths_um=mo["depths_um"]),
            shanks=[dict(shank=s["shank"], dead=s["dead"], uv_per_count=s["uv_per_count"]) for s in shanks]))
        print(f"probe {p}: {n} samples, tiles {t_first} to {t_last}")
    man = dict(
        version=2,
        description="Raw voltage for the raster viewer's underlay, FD_016 2026-11-13/1; see "
                    "raw_underlay/make_raw.py in github.com/SteinmetzLab/rasterViewer",
        fs=FS, tile=TILE, margin=MARGIN, n_channels=384, row_pitch_um=15.0,
        highpass=dict(hz=HP_HZ, order=HP_ORDER, sos=sos.tolist()),
        variant=dict(key=f"bits{SHIFT}", shift=SHIFT, scale=float(2 ** SHIFT), offset=(2 ** SHIFT - 1) / 2,
                     label=f"{SHIFT} lowest bits dropped"),
        path="p{probe}/sh{shank}/0.{tile}", probes=probes)
    json.dump(man, open(out / "manifest.json", "w"), indent=1)
    json.dump(checks, open(out / "clock_check.json", "w"), indent=1)
    print("manifest written")


# ---------------------------------------------------------------------------- convert
def s3_client(keys_path):
    import boto3
    from botocore.config import Config
    k = json.load(open(keys_path))
    cl = boto3.client("s3", endpoint_url=f"https://{k['account_id'].strip()}.r2.cloudflarestorage.com",
                      aws_access_key_id=k["access_key_id"].strip(),
                      aws_secret_access_key=k["secret_access_key"].strip(), region_name="auto",
                      config=Config(retries={"max_attempts": 8, "mode": "adaptive"}, max_pool_connections=16))
    return cl, k["bucket"]


def encode_tile(v):
    """Stored bytes of one tile: v is int16 (384, TILE); delta over the flattened array, zlib."""
    f = v.reshape(-1)
    d = f.copy()
    d[1:] = f[1:] - f[:-1]                                  # int16 arithmetic wraps, as intended
    return zlib.compress(d.astype("<i2").tobytes(), ZLIB_LEVEL)


def decode_tile(buf):
    """Inverse of encode_tile -> int16 (384, TILE)."""
    d = np.frombuffer(zlib.decompress(buf), "<i2")
    assert d.size == 384 * TILE
    return np.cumsum(d, dtype=np.int16).reshape(384, TILE)


_W = {}


def _init(session, out, keys):
    _W.update(session=session, out=Path(out), cb={}, s3=None, pool=None)
    if keys:
        _W["s3"], _W["bucket"] = s3_client(keys)
        _W["pool"] = ThreadPoolExecutor(5)


def _put(key, body):
    try:
        _W["s3"].put_object(Bucket=_W["bucket"], Key=key, Body=body,
                            ContentType="application/octet-stream", CacheControl=CACHE_CHUNK)
        return True
    except Exception as e:                                   # `sync` uploads what is missing
        sys.stderr.write(f"upload of {key} failed: {type(e).__name__}: {e}\n")
        return False


def _chunk(task):
    """One second of one shank's .cbin -> its tiles (those in `tiles`), written and uploaded."""
    p, sh, c, tiles = task
    cb = _W["cb"].get((p, sh))
    if cb is None:
        cb = _W["cb"][(p, sh)] = Cbin(raw_path(_W["session"], p, sh))
    blk = cb.chunk(c)[:, :384]                               # samples x sites
    res, puts = [], []
    for T in tiles:
        a = T * TILE - c * FS
        seg = blk[a:a + TILE]
        if seg.shape[0] < TILE:                              # the recording's last, partial tile
            seg = np.concatenate([seg, np.zeros((TILE - seg.shape[0], 384), np.int16)])
        v = np.ascontiguousarray((seg >> SHIFT).T)
        buf = encode_tile(v)
        rel = f"p{p}/sh{sh}/0.{T}"
        (_W["out"] / rel).write_bytes(buf)
        res.append(dict(T=int(T), n=len(buf), md5=hashlib.md5(buf).hexdigest(),
                        sum=int(v.sum(dtype=np.int64)), up=False))
        if _W["s3"]:
            puts.append(_W["pool"].submit(_put, PREFIX + rel, buf))
    for r, f in zip(res, puts):
        r["up"] = f.result()
    return p, sh, c, res


def read_progress(out, p, sh):
    f = Path(out) / "progress" / f"p{p}_sh{sh}.jsonl"
    done = {}
    if f.exists():
        for line in open(f):
            line = line.strip()
            if line:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:                 # a line cut short by a crash
                    continue
                done[r["T"]] = r
    return done


def tile_range(man, pr, t0=None, t1=None):
    """First and last tile to write for a probe; with t0 / t1 (viewer s), those around them."""
    lo, hi = pr["tile_first"], pr["tile_last"]
    if t0 is not None:
        k = np.fromfile(Path(man["_dir"]) / pr["clock"], "<f8").reshape(-1, 2)
        lo = max(lo, int(np.interp(t0, k[:, 1], k[:, 0]) // TILE) - 1)
        hi = min(hi, int(np.interp(t1, k[:, 1], k[:, 0]) // TILE) + 1)
    return lo, hi


def do_convert(a):
    import multiprocessing as mp
    out = Path(a.out)
    man = json.load(open(out / "manifest.json"))
    man["_dir"] = str(out)
    (out / "progress").mkdir(exist_ok=True)
    t_begin = time.time()
    n_done = n_bytes = n_fail = 0
    with mp.Pool(a.workers, initializer=_init, initargs=(a.session, str(out), a.keys)) as pool:
        for pr in man["probes"]:
            p = pr["probe"]
            if a.probes is not None and p not in a.probes:
                continue
            lo, hi = tile_range(man, pr, a.t0, a.t1)
            per = FS // TILE
            for s in pr["shanks"]:
                sh = s["shank"]
                done = read_progress(out, p, sh)
                tasks = []
                for c in range(lo // per, hi // per + 1):
                    tiles = [T for T in range(c * per, (c + 1) * per) if lo <= T <= hi and T not in done]
                    if tiles:
                        tasks.append((p, sh, c, tiles))
                if not tasks:
                    print(f"probe {p} shank {sh}: nothing to do", flush=True)
                    continue
                t_sh = time.time()
                b_sh = 0
                with open(out / "progress" / f"p{p}_sh{sh}.jsonl", "a") as log:
                    for k, (_, _, c, res) in enumerate(pool.imap_unordered(_chunk, tasks)):
                        for r in res:
                            log.write(json.dumps(r) + "\n")
                            n_done += 1
                            b_sh += r["n"]
                            n_fail += (a.keys is not None and not r["up"])
                        if k % 50 == 0:
                            log.flush()
                n_bytes += b_sh
                nt = sum(len(t[3]) for t in tasks)
                print(f"probe {p} shank {sh}: {nt} tiles in {time.time() - t_sh:.0f} s, "
                      f"{nt * 384 * TILE * 2 / b_sh:.2f}x smaller than int16; total {n_done} tiles, "
                      f"{n_bytes / 1e9:.1f} GB, {n_fail} failed uploads, {(time.time() - t_begin) / 60:.1f} min",
                      flush=True)
    print(f"convert finished: {n_done} tiles, {n_bytes / 1e9:.2f} GB, {n_fail} failed uploads, "
          f"{(time.time() - t_begin) / 60:.1f} min", flush=True)


# ---------------------------------------------------------------------------- refs
def do_refs(a):
    """Reference tiles for the viewer's ?rawtest=1, one chunk sum per shank, and the sizes."""
    out = Path(a.out)
    man = json.load(open(out / "manifest.json"))
    (out / "test").mkdir(exist_ok=True)
    v = man["variant"]
    tests, n_tiles, n_bytes = [], 0, 0
    for pr in man["probes"]:
        p = pr["probe"]
        k = np.fromfile(out / pr["clock"], "<f8").reshape(-1, 2)
        T = int(np.interp(a.t_ref, k[:, 1], k[:, 0]) // TILE)
        for s in pr["shanks"]:
            done = read_progress(out, p, s["shank"])
            n_tiles += len(done)
            n_bytes += sum(r["n"] for r in done.values())
            if T in done:
                s["check"] = dict(tile=T, sum=done[T]["sum"])
        if p in a.ref_probes:
            sh = a.ref_shank
            raw = np.concatenate([decode_tile((out / f"p{p}/sh{sh}/0.{t}").read_bytes()) for t in (T - 1, T, T + 1)], axis=1)
            s = pr["shanks"][sh]
            ref = process_tile(raw, 1, s["dead"], s["uv_per_count"], v["scale"], v["offset"])
            f = f"test/ref_p{p}_sh{sh}_tile{T}.f32"
            ref.tofile(out / f)
            tests.append(dict(probe=p, shank=sh, tile=T, file=f))
    man["test"] = tests
    if n_tiles:
        v.update(tiles_written=n_tiles, stored_bytes=n_bytes, ratio=n_tiles * 384 * TILE * 2 / n_bytes)
    json.dump(man, open(out / "manifest.json", "w"), indent=1)
    print(f"{len(tests)} reference tiles; {n_tiles} tiles, {n_bytes / 1e9:.2f} GB, ratio {v.get('ratio')}")


# ---------------------------------------------------------------------------- sync
def do_sync(a):
    out = Path(a.out)
    man = json.load(open(out / "manifest.json"))
    s3, bucket = s3_client(a.keys)
    have = {}
    for page in s3.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=PREFIX):
        for o in page.get("Contents", []):
            have[o["Key"]] = (o["Size"], o["ETag"].strip('"'))
    print(f"{len(have)} objects in the bucket under {PREFIX}")
    want, expected = [], set()
    for pr in man["probes"]:
        for s in pr["shanks"]:
            rel = f"p{pr['probe']}/sh{s['shank']}"
            for T, r in read_progress(out, pr["probe"], s["shank"]).items():
                key = f"{PREFIX}{rel}/0.{T}"
                expected.add(key)
                if have.get(key) != (r["n"], r["md5"]):
                    want.append((key, out / rel / f"0.{T}", r["md5"]))
    print(f"{len(expected)} tiles in the progress files, {len(want)} to upload")

    def up(item):
        key, path, md5 = item
        buf = path.read_bytes()
        if hashlib.md5(buf).hexdigest() != md5:
            return f"{path}: file does not match the progress file"
        s3.put_object(Bucket=bucket, Key=key, Body=buf, ContentType="application/octet-stream",
                      CacheControl=CACHE_CHUNK)
        return None
    with ThreadPoolExecutor(16) as ex:
        bad = [e for e in ex.map(up, want) if e]
    for e in bad[:20]:
        print("  ", e)
    small = ["manifest.json"] + [pr["clock"] for pr in man["probes"]] + [pr["motion"]["file"] for pr in man["probes"]] \
        + [f"p{pr['probe']}/sh{s['shank']}/.zarray" for pr in man["probes"] for s in pr["shanks"]] \
        + [t["file"] for t in man.get("test", [])]
    for rel in small:
        ct = "application/json" if rel.endswith((".json", ".zarray")) else "application/octet-stream"
        s3.put_object(Bucket=bucket, Key=PREFIX + rel, Body=(out / rel).read_bytes(), ContentType=ct,
                      CacheControl=CACHE_SMALL)
        expected.add(PREFIX + rel)
    stray = sorted(set(have) - expected)
    print(f"uploaded {len(want) - len(bad)} tiles ({len(bad)} problems) and {len(small)} small files; "
          f"{len(stray)} objects in the bucket are not expected{': ' + ', '.join(stray[:5]) if stray else ''}")
    return len(bad)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["meta", "convert", "refs", "sync"])
    ap.add_argument("--session", default="Y:/Subjects/FD_016/2026-11-13/1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--keys", default=None, help="JSON with the bucket's credentials; convert uploads if given")
    ap.add_argument("--probes", type=int, nargs="*", default=None)
    ap.add_argument("--t0", type=float, default=None, help="convert only around this stretch (viewer s)")
    ap.add_argument("--t1", type=float, default=None)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--t-ref", type=float, default=1125.05, help="viewer time of the test tiles")
    ap.add_argument("--ref-probes", type=int, nargs="*", default=[0, 1, 4])
    ap.add_argument("--ref-shank", type=int, default=0)
    a = ap.parse_args()
    assert (a.t0 is None) == (a.t1 is None)
    r = {"meta": do_meta, "convert": do_convert, "refs": do_refs, "sync": do_sync}[a.step](a)
    sys.exit(1 if r else 0)
