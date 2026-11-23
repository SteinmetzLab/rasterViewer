"""Does the underlay put each sorted spike on its waveform? (time and depth)

    python validate_alignment.py --raw Y:/temp/fd016_raw_underlay [--shank 0]

Uses only what the viewer uses: the lossless Zarr array, the clock pairs (clock_p<probe>.f64)
and the viewer's spike table (ALF times, unit depths from data/meta.json). For every QC unit of
the shank with at least 30 spikes in the stretch:

  time   The viewer time of each spike is mapped to a probe sample through the clock pairs
         (linear between pairs, as raw_underlay.js does). The spike-triggered average of the
         high-passed voltage is taken around that sample on every site within 150 um of the
         unit. Lag = sample of the deepest trough on the largest site, relative to the mapped
         sample. For comparison the same is done with a single straight line fitted to the
         clock pairs, which is what a map without the sync wobble would give.
  depth  Center of the trough amplitude over sites (weights = amplitude squared) against the
         unit's depth in the viewer.

Writes validate_alignment.json next to the raw data and prints a summary.
"""
import argparse
import json
from pathlib import Path

import numpy as np

FS = 30000
HALF = 45        # samples on each side of the spike in the average


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--shank", type=int, default=0)
    ap.add_argument("--data", default=str(Path(__file__).resolve().parent.parent / "data"))
    ap.add_argument("--alf", default="D:/temp/fd016_share/alf")
    a = ap.parse_args()
    import zarr
    from scipy import signal

    raw = Path(a.raw)
    man = json.load(open(raw / "manifest.json"))
    probe = man["probe"]
    sh = next(s for s in man["shanks"] if s["shank"] == a.shank)
    pairs = np.fromfile(raw / man["clock"], "<f8").reshape(-1, 2)
    cs, ct = pairs[:, 0], pairs[:, 1]
    lin = np.polyfit(ct, cs, 1)                       # straight-line clock, for comparison
    z = zarr.open_array(str(raw / f"p{probe}" / f"sh{a.shank}" / "lossless"), mode="r")
    x = z[:]                                          # channels x samples, int16
    n = x.shape[1]
    sos = np.asarray(man["highpass"]["sos"])
    y_site = np.asarray(sh["y_um"])
    uv = sh["uv_per_count"]

    meta = json.load(open(Path(a.data) / "meta.json"))
    rows = meta["rows"]
    st = np.load(Path(a.alf) / "spikes.times.npy")
    sc = np.load(Path(a.alf) / "spikes.clusters.npy")
    t0, t1 = np.interp([man["s0"] + HALF + 1, man["s0"] + n - HALF - 2], cs, ct)
    i0, i1 = np.searchsorted(st, [t0, t1])
    st, sc = st[i0:i1], sc[i0:i1]

    units = [i for i in range(meta["n_rows"]) if rows["probe"][i] == probe and rows["shank"][i] == a.shank]
    spikes = {}
    for i in units:
        t = st[sc == rows["cluster"][i]]
        if t.size >= 30:
            spikes[i] = t
    print(f"shank {a.shank}: {len(units)} QC units, {len(spikes)} with 30 or more spikes in the stretch")

    # spike-triggered averages, one site at a time (each site is filtered once)
    sta = {i: {} for i in spikes}
    sta_lin = {i: {} for i in spikes}
    idx = np.arange(-HALF, HALF + 1)
    for c in range(x.shape[0]):
        near = [i for i in spikes if abs(rows["y_um"][i] - y_site[c]) <= 150]
        if not near or c in sh["dead"]:
            continue
        v = signal.sosfiltfilt(sos, x[c].astype(np.float64)) * uv
        for i in near:
            s = np.rint(np.interp(spikes[i], ct, cs) - man["s0"]).astype(int)
            sta[i][c] = v[s[:, None] + idx[None, :]].mean(0)
            s2 = np.clip(np.rint(np.polyval(lin, spikes[i]) - man["s0"]).astype(int), HALF, n - HALF - 1)
            sta_lin[i][c] = v[s2[:, None] + idx[None, :]].mean(0)

    out = []
    for i in spikes:
        chans = sorted(sta[i])
        W = np.array([sta[i][c] for c in chans])               # sites x samples
        W2 = np.array([sta_lin[i][c] for c in chans])
        amp = -W.min(axis=1)
        k = int(np.argmax(amp))
        lag = int(np.argmin(W[k])) - HALF
        w = np.maximum(amp, 0) ** 2
        yc = float((w * y_site[chans]).sum() / w.sum())
        out.append(dict(row=i, cluster=rows["cluster"][i], n_spikes=int(spikes[i].size),
                        amp_uv=float(amp[k]), lag_samples=lag,
                        amp_straight_line_uv=float(-W2[k].min()),
                        y_viewer_um=float(rows["y_um"][i]), y_raw_um=yc,
                        y_peak_site_um=float(y_site[chans[k]])))
    big = [o for o in out if o["amp_uv"] >= 30]
    lag = np.array([o["lag_samples"] for o in big])
    dy = np.array([o["y_raw_um"] - o["y_viewer_um"] for o in big])
    ratio = np.array([o["amp_straight_line_uv"] / o["amp_uv"] for o in big])
    summ = dict(shank=a.shank, n_units=len(out), n_units_30uv=len(big),
                lag_samples=dict(median=float(np.median(lag)), p5=float(np.percentile(lag, 5)),
                                 p95=float(np.percentile(lag, 95)), within_2=float(np.mean(np.abs(lag) <= 2))),
                depth_um=dict(median=float(np.median(dy)), median_abs=float(np.median(np.abs(dy))),
                              p95_abs=float(np.percentile(np.abs(dy), 95))),
                straight_line_amp_ratio=dict(median=float(np.median(ratio)), p5=float(np.percentile(ratio, 5))))
    json.dump(dict(summary=summ, units=out), open(raw / f"validate_alignment_sh{a.shank}.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
