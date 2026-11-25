"""Check tiles of the display copy (make_display.py) against the original recording, with code
that shares nothing with the converter: the recording is read with the mtscomp package, the
filter is scipy's sosfiltfilt over a longer stretch, and local tiles are read through
zarr-python.

    python verify_display.py --work D:/temp/fd016_disp [--url https://.../d1/] [--n 3]
                             [--probes 0] [--local]

For every shank, n tiles of each level at random among those the progress files list (plus the
first and last of the range when they are there). With --url the tiles are downloaded from the
bucket (only tiles recorded as uploaded are picked); with --local they are read from the work
folder through zarr-python (only tiles whose file is there are picked: a missing chunk would
read as zeros).

What may differ. The converter filters each second with 10 ms of margin, this script filters
5 s around the tile in one piece, so a value can land on the other side of a rounding boundary:
a difference of one count is allowed where the unrounded value is within TOL of the boundary,
and counted. Anything else is an error. Level 1 is also checked against level 0 recomputed
here: each value must be the largest-magnitude sample of its block.
"""
import argparse
import json
import subprocess
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy import signal

FS = 30000
TOL = 0.02          # of a count: how far from a rounding boundary a one-count difference may be


def fetch(url):
    r = subprocess.run(["curl", "-s", "-f", url], capture_output=True)     # curl: Python's TLS store is stale here
    if r.returncode:
        raise RuntimeError(f"{url}: curl exit {r.returncode}")
    return r.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--session", default="Y:/Subjects/FD_016/2026-11-13/1")
    ap.add_argument("--url", default=None)
    ap.add_argument("--local", action="store_true")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--probes", type=int, nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    assert a.url or a.local
    import mtscomp
    work = Path(a.work)
    man = json.load(open(work / "manifest.json"))
    raw = json.load(open(work / "v1" / "manifest.json"))
    step, clip = man["step_uv"], man["clip"]
    sos = signal.butter(man["highpass"]["order"], man["highpass"]["hz"], "hp", fs=FS, output="sos")
    rng = np.random.default_rng(a.seed)
    n_tiles = n_values = n_flip = n_bad = n_none = 0
    for pr in man["probes"]:
        p = pr["probe"]
        if a.probes is not None and p not in a.probes:
            continue
        run, imec = pr["name"].split("_")
        for s in pr["shanks"]:
            sh = s["shank"]
            rs = raw["probes"][p]["shanks"][sh]
            done = {}
            f = work / "progress" / f"p{p}_sh{sh}.jsonl"
            for line in open(f) if f.exists() else []:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                for t in r["tiles"]:
                    done[(t["L"], t["T"])] = t
            cb = mtscomp.decompress(Path(a.session) / "raw_ephys_data" / f"{run}_g0_{imec}" /
                                    f"{run}_g0_t0.{imec}.sh{sh}.ap.cbin")
            good = np.setdiff1d(np.arange(384), rs["dead"])
            for lv, rg in zip(man["levels"], pr["levels"]):
                L, k, nt = lv["level"], lv["block"], lv["tile"]
                tdir = work / f"p{p}/sh{sh}/L{L}"
                have = sorted(T for (l, T), t in done.items() if l == L and
                              (t["up"] if a.url else (tdir / f"0.{T}").exists()))
                if not have:
                    print(f"probe {p} shank {sh} level {L}: no tile to check")
                    n_none += 1
                    continue
                pick = set(rng.choice(have, min(a.n, len(have)), replace=False).tolist())
                pick |= {T for T in (rg["tile_first"], rg["tile_last"]) if T in have}
                for T in sorted(pick):
                    rel = lv["path"].replace("{probe}", str(p)).replace("{shank}", str(sh)).replace("{tile}", str(T))
                    if a.url:
                        buf = fetch(a.url + rel)
                        got = np.frombuffer(zlib.decompress(buf), np.int8).reshape(384, nt)
                        if len(buf) != done[(L, T)]["n"]:
                            print(f"probe {p} shank {sh} level {L} tile {T}: size differs from the progress file")
                            n_bad += 1
                    else:
                        import zarr
                        z = zarr.open(str(work / f"p{p}/sh{sh}/L{L}"), mode="r")
                        got = np.asarray(z[:, T * nt:(T + 1) * nt])
                    if int(got.sum(dtype=np.int64)) != done[(L, T)]["sum"]:
                        print(f"probe {p} shank {sh} level {L} tile {T}: sum differs from the progress file")
                        n_bad += 1
                    s0, s1 = T * nt * k, (T + 1) * nt * k
                    a0, a1 = max(0, s0 - 2 * FS), min(cb.shape[0], s1 + 2 * FS)
                    x = np.asarray(cb[a0:a1, :384], np.float64).T
                    x -= x.mean(axis=1, keepdims=True)
                    y = signal.sosfiltfilt(sos, x, axis=1, padtype=None)
                    y[rs["dead"]] = 0.0
                    y[good] -= np.median(y[good], axis=0, keepdims=True)
                    u = np.zeros((384, s1 - s0))
                    u[:, :min(s1, a1) - s0] = y[:, s0 - a0:s1 - a0] * rs["uv_per_count"] / step
                    if k > 1:
                        b = u.reshape(384, nt, k)
                        mag = np.abs(b).max(axis=2)
                        # the stored value must be a sample of its block (within rounding) and
                        # as large in magnitude as the largest
                        near = np.abs(b - got[..., None]).min(axis=2)
                        want_mag = np.clip(np.rint(mag), 0, clip)
                        d = np.abs(np.abs(got.astype(np.int64)) - want_mag)
                        edge = np.abs(mag - np.floor(mag) - 0.5) < TOL
                        bad = (d > 1) | ((d == 1) & ~edge) | ((near > 0.5 + TOL) & (np.abs(got) < clip))
                        flip = (d == 1) & edge
                    else:
                        want = np.clip(np.rint(u), -clip, clip)
                        d = np.abs(got - want)
                        edge = np.abs(u - np.floor(u) - 0.5) < TOL
                        bad = (d > 1) | ((d == 1) & ~edge)
                        flip = (d == 1) & edge
                    n_tiles += 1; n_values += got.size; n_flip += int(flip.sum()); n_bad += int(bad.sum())
                    if bad.any():
                        print(f"probe {p} shank {sh} level {L} tile {T}: {int(bad.sum())} values differ")
            print(f"probe {p} shank {sh}: checked", flush=True)
    print(f"{n_tiles} tiles, {n_values} values: {n_bad} differ; {n_flip} differ by one count within "
          f"{TOL} of a rounding boundary (allowed); {n_none} shank-levels had no tile to check")
    return 1 if n_bad or n_none or not n_tiles else 0


if __name__ == "__main__":
    sys.exit(main())
