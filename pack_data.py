"""Compress the spike files for the web (about 366 MB -> 115 MB), and read them back.

    python pack_data.py RAW_DIR      # RAW_DIR holds ticks.u32 and rows.u16 (from prep_data.py)

Writes into data/ next to this file:
  spikes_ticks.gz.NN  spike times: the first tick, then each spike's tick minus the previous
                      one as uint8 (255 = escape: the true delta is the next value in
                      manifest "tick_escapes"), gzip level 9, split into parts
  spikes_rows.gz.NN   each spike's display row as uint16, bytes split into two planes (all
                      high bytes, then all low bytes), gzip level 9, split into parts
  manifest.json       part names and sizes, spike count, escapes, and checksums of the decoded
                      arrays, which the viewer verifies after loading

Concatenating a stream's parts in order gives one gzip file. Parts stay under 45 MB, below
GitHub's 50 MB warning and 100 MB limit. Deltas are small because spikes are time-sorted
(61 M spikes in 73 M ticks; the largest delta in FD_016 is 35 ticks).

load_spikes(data_dir) reads either the raw files or the packed ones, for reference.py.
"""
import gzip
import json
import sys
import zlib
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PART = 45_000_000


def checksums(ticks, rows):
    """The same numbers viewer.js computes after decoding (exact in float64)."""
    t = ticks.astype(np.int64); r = rows.astype(np.int64)
    i = np.arange(len(r), dtype=np.int64)
    return dict(sum_ticks=int(t.sum()), sum_rows=int(r.sum()),
                sum_rows_x_index_mod_1024=int((r * (i & 1023)).sum()),
                first_tick=int(t[0]), last_tick=int(t[-1]))


def pack(raw_dir, out_dir=HERE / "data"):
    ticks = np.fromfile(Path(raw_dir) / "ticks.u32", "<u4")
    rows = np.fromfile(Path(raw_dir) / "rows.u16", "<u2")
    assert len(ticks) == len(rows) and np.all(np.diff(ticks.astype(np.int64)) >= 0)
    d = np.diff(ticks.astype(np.int64))
    esc = d >= 255
    d8 = np.where(esc, 255, d).astype(np.uint8)
    planes = rows.view(np.uint8).reshape(-1, 2)
    streams = {"ticks": d8.tobytes(),
               "rows": np.ascontiguousarray(planes[:, 1]).tobytes()
                       + np.ascontiguousarray(planes[:, 0]).tobytes()}
    for f in out_dir.glob("spikes_*.gz.*"):
        f.unlink()
    man = dict(format="rasterViewer-spikes-1", n_spikes=int(len(ticks)),
               tick0=int(ticks[0]), tick_escapes=[int(x) for x in d[esc]],
               checksums=checksums(ticks, rows), streams={})
    for name, raw in streams.items():
        comp = gzip.compress(raw, compresslevel=9, mtime=0)
        parts = []
        for k in range(0, len(comp), PART):
            fn = f"spikes_{name}.gz.{k // PART:02d}"
            (out_dir / fn).write_bytes(comp[k:k + PART])
            parts.append(dict(file=fn, bytes=len(comp[k:k + PART])))
        man["streams"][name] = dict(raw_bytes=len(raw), parts=parts)
        print(f"{name}: {len(raw) / 1e6:.0f} MB -> {len(comp) / 1e6:.1f} MB in {len(parts)} parts")
    json.dump(man, open(out_dir / "manifest.json", "w"), indent=1)
    t2, r2 = load_spikes(out_dir)          # round trip
    assert np.array_equal(t2, ticks) and np.array_equal(r2, rows), "round trip failed"
    print("round trip exact;", man["checksums"])


def load_spikes(data_dir=HERE / "data"):
    """(ticks uint32, rows uint16), from raw files if present, else from the packed parts."""
    data_dir = Path(data_dir)
    if (data_dir / "ticks.u32").exists():
        return (np.fromfile(data_dir / "ticks.u32", "<u4"),
                np.fromfile(data_dir / "rows.u16", "<u2"))
    man = json.load(open(data_dir / "manifest.json"))
    raw = {}
    for name, s in man["streams"].items():
        comp = b"".join((data_dir / p["file"]).read_bytes() for p in s["parts"])
        raw[name] = zlib.decompress(comp, 16 + zlib.MAX_WBITS)
        assert len(raw[name]) == s["raw_bytes"]
    n = man["n_spikes"]
    d = np.frombuffer(raw["ticks"], np.uint8).astype(np.int64)
    esc = d == 255
    d[esc] = man["tick_escapes"]
    ticks = np.empty(n, np.int64); ticks[0] = man["tick0"]; np.cumsum(d, out=ticks[1:])
    ticks[1:] += man["tick0"]
    b = np.frombuffer(raw["rows"], np.uint8)
    rows = (b[:n].astype(np.uint16) << 8) | b[n:]
    ticks = ticks.astype(np.uint32); rows = rows.astype(np.uint16)
    assert checksums(ticks, rows) == man["checksums"], "checksum mismatch"
    return ticks, rows


if __name__ == "__main__":
    pack(sys.argv[1] if len(sys.argv) > 1 else HERE / "data")
