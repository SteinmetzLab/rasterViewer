"""Write the viewer's data files for FD_016 2026-11-13/1 (run once, ~1 min).

Inputs:  D:/temp/fd016_share/alf (spikes.times, spikes.clusters; 4,325 QC units)
         D:/temp/fd016_raster/units.npz (meta: probe, shank, median Kilosort spike y per unit,
         from the Breaking Barriers movie export; its unit index equals the ALF cluster index,
         checked 2026-11-21: per-unit counts in 1102-1182 s agree exactly)
Outputs: data/ next to this file (then run pack_data.py data to make the packed copies the
         viewer and the repository use)
  ticks.u32  spike times as 30 kHz sample counts on the Timeline clock, round(t * 30000),
             sorted ascending (uint32, little-endian)
  rows.u16   display row of each spike (uint16), same order
  meta.json  per-row table (cluster, probe, shank, depth, session mean rate, shank segment
             bounds), time extent, colormap tables, 1 s population rate for the overview

Display rows: probe, then shank, then y descending (shallow at top), the movie's order. y is the
median Kilosort spike position along the shank, in um above the deepest site (larger = closer
to the brain surface); it is not depth below the surface, which needs histology.
Each (probe, shank) is one "segment"; the viewer never averages across segment boundaries.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib

ALF = Path(r"D:/temp/fd016_share/alf")
MOVIE_NPZ = Path(r"D:/temp/fd016_raster/units.npz")
OUT = Path(__file__).resolve().parent / "data"
FS = 30000
T_START, T_END = 1.9685, 2430.811     # synced period (ALF README), seconds
PROBE_HUES = ["#14b899", "#ff6f4e", "#4ea8ff", "#f2c14e", "#b07cff", "#a4d65e", "#ff5fa2"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    st = np.load(ALF / "spikes.times.npy")
    sc = np.load(ALF / "spikes.clusters.npy")
    meta = np.load(MOVIE_NPZ)["meta"]
    n = meta.shape[0]
    assert sc.max() == n - 1
    assert np.all(np.diff(st) >= 0), "spike times must be sorted"
    order = np.lexsort((-meta[:, 2], meta[:, 1], meta[:, 0]))   # row -> cluster
    row_of = np.empty(n, np.int64); row_of[order] = np.arange(n)

    ticks = np.round(st * FS).astype(np.uint32)
    assert np.all(np.diff(ticks.astype(np.int64)) >= 0)
    rows = row_of[sc].astype(np.uint16)
    ticks.tofile(OUT / "ticks.u32")
    rows.tofile(OUT / "rows.u16")

    probe = meta[order, 0].astype(int); shank = meta[order, 1].astype(int)
    seg_id = probe * 4 + shank
    assert np.all(np.diff(seg_id) >= 0)
    seg_start = np.zeros(n, int); seg_end = np.zeros(n, int)
    segments = []
    for s in np.unique(seg_id):
        idx = np.flatnonzero(seg_id == s)
        seg_start[idx] = idx[0]; seg_end[idx] = idx[-1] + 1
        segments.append(dict(probe=int(s // 4), shank=int(s % 4), start=int(idx[0]),
                             end=int(idx[-1] + 1)))
    dur = T_END - T_START
    r0 = np.bincount(rows, minlength=n) / dur

    pop = np.bincount(np.floor(st).astype(int), minlength=int(np.ceil(T_END)) + 1)
    cmaps = {}
    for name in ["magma", "inferno", "viridis", "gray", "RdBu_r", "PuOr_r"]:
        lut = (matplotlib.colormaps[name](np.linspace(0, 1, 256))[:, :3] * 255).round()
        cmaps[name] = lut.astype(int).ravel().tolist()

    out = dict(
        dataset="FD_016 2026-11-13/1, 4,325 QC units on 7 NP2.0 Quad-Base probes",
        fs=FS, n_rows=n, n_spikes=int(len(ticks)), t_start=T_START, t_end=T_END,
        rows=dict(cluster=order.tolist(), probe=probe.tolist(), shank=shank.tolist(),
                  y_um=np.round(meta[order, 2], 1).tolist(),   # height above the shank's deepest site
                  r0=np.round(r0, 6).tolist(),
                  seg_start=seg_start.tolist(), seg_end=seg_end.tolist()),
        segments=segments, probe_hues=PROBE_HUES,
        pop_rate_1s=pop.tolist(), colormaps=cmaps,
    )
    json.dump(out, open(OUT / "meta.json", "w"))
    print(f"{len(ticks):,} spikes, {n} rows, {len(segments)} segments; "
          f"ticks {ticks[0]}..{ticks[-1]}; written to {OUT}")


if __name__ == "__main__":
    main()
