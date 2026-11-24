"""Independent check of the arrays written by make_raw.py.

    python verify_raw.py --out Y:/temp/fd016_raw_v1 [--n 3] [--url https://...r2.dev/v1/]

For --n random tiles of every shank (among those listed in the progress files):
  - the tile is read through zarr-python, i.e. by a standard Zarr reader and not by the
    converter's own code (with --url, also from the bucket's public address);
  - the same samples are read from the lab's .cbin file with the mtscomp package (not the
    converter's cbin.py);
  - stored value == sample >> shift for every sample, and the value the viewer uses
    (value * scale + offset) is within (2**shift - 1) / 2 counts of the sample.
Tiles around the recording's ends are always included. Needs zarr, numcodecs, mtscomp.
"""
import argparse
import json
import sys
import urllib.request
import zlib
from pathlib import Path

import numpy as np
import zarr
from mtscomp import Reader

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_raw import read_progress                      # noqa: E402  (the list of tiles only)
from make_raw_proto import raw_path                     # noqa: E402


def fetch_tile(url):
    req = urllib.request.Request(url, headers={"User-Agent": "rasterViewer-verify"})
    import ssl
    import certifi
    with urllib.request.urlopen(req, timeout=60, context=ssl.create_default_context(cafile=certifi.where())) as f:
        d = np.frombuffer(zlib.decompress(f.read()), "<i2")
    return np.cumsum(d, dtype=np.int16).reshape(384, -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--session", default="Y:/Subjects/FD_016/2026-11-13/1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--url", default=None, help="also compare with the hosted copy at this address")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    man = json.load(open(out / "manifest.json"))
    tile, v = man["tile"], man["variant"]
    rng = np.random.default_rng(a.seed)
    n_tiles = n_bad = n_url = 0
    max_err = 0.0
    for pr in man["probes"]:
        for s in pr["shanks"]:
            p, sh = pr["probe"], s["shank"]
            done = sorted(read_progress(out, p, sh))
            if not done:
                continue
            pick = sorted(set([done[0], done[-1]] + rng.choice(done, min(a.n, len(done)), replace=False).tolist()))
            z = zarr.open_array(str(out / f"p{p}" / f"sh{sh}"), mode="r")
            assert z.shape == (384, pr["n_samples"]) and z.chunks == (384, tile) and z.dtype == np.int16
            cb = raw_path(a.session, p, sh)
            r = Reader()
            r.open(cb, cb.with_suffix(".ch"))
            for T in pick:
                a0, a1 = T * tile, min((T + 1) * tile, pr["n_samples"])
                got = np.asarray(z[:, a0:a1])
                orig = np.asarray(r[a0:a1])[:, :384].T
                ok = np.array_equal(got, orig >> v["shift"])
                err = np.abs(got.astype(np.float64) * v["scale"] + v["offset"] - orig).max()
                max_err = max(max_err, err)
                if a.url:
                    web = fetch_tile(a.url.rstrip("/") + "/" + man["path"].replace("{probe}", str(p))
                                     .replace("{shank}", str(sh)).replace("{tile}", str(T)))
                    ok = ok and np.array_equal(web[:, :a1 - a0], got)
                    n_url += 1
                n_tiles += 1
                n_bad += not ok
                if not ok:
                    print(f"MISMATCH probe {p} shank {sh} tile {T}")
            r.close()
            print(f"probe {p} shank {sh}: tiles {pick} checked", flush=True)
    print(f"{n_tiles} tiles checked ({n_url} also from the bucket), {n_bad} differ; largest error of the "
          f"viewer's value {max_err} counts (limit {(2 ** v['shift'] - 1) / 2})")
    sys.exit(1 if n_bad or max_err > (2 ** v["shift"] - 1) / 2 else 0)


if __name__ == "__main__":
    main()
