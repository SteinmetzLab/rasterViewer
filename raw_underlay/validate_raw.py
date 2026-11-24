"""Does the underlay put each sorted spike on its waveform? Time, depth, and the motion shift.

    python validate_raw.py --raw Y:/temp/fd016_raw_v1 --probe 1 --shank 0 --t0 5 --t1 70

Uses only what the viewer uses: the stored tiles (value * scale + offset), the clock knots
(clock_p<probe>.f64), the motion estimate (motion_p<probe>.f32) and the viewer's spike table
(ALF times, unit depths from data/meta.json). For every QC unit of the shank with at least 30
spikes in the stretch:

  time   The viewer time of each spike is mapped to a probe sample through the clock knots
         (linear between knots, as raw_underlay.js does). The spike-triggered average of the
         high-passed voltage is taken around that sample on every site within 150 um of the
         unit. Lag = sample of the deepest trough on the largest site, relative to the mapped
         sample. For comparison the same is done with a single straight line fitted to the
         knots, which is what a map without the sync wobble would give.
  depth  Center of the trough amplitude over sites (weights = amplitude squared), y_raw, in
         the probe's own coordinates. The unit's depth in the viewer, y_unit, is in the sort's
         motion-corrected coordinates. d = the motion estimate at the unit's depth, averaged
         over its spikes (per spike: the 1 s bin holding its sample, linear between the depth
         windows, constant beyond the outer ones). Three candidates for y_unit are compared:
         y_raw - d (the viewer's rule: a site at y is drawn at y - d), y_raw, and y_raw + d.

Writes validate_raw_p<probe>_sh<shank>.json next to the raw data and prints a summary.
"""
import argparse
import json
import zlib
from pathlib import Path

import numpy as np

HALF = 45        # samples on each side of the spike in the average


def displacement(mo, disp, sample, y, fs):
    """Motion estimate (um) at probe samples and depths, by the viewer's rule."""
    b = np.clip(np.floor((sample / fs - mo["t0_s"]) / mo["dt_s"] + 0.5).astype(int), 0, disp.shape[0] - 1)
    dep = np.asarray(mo["depths_um"])
    u = np.clip((y - dep[0]) / (dep[1] - dep[0]), 0, dep.size - 1)
    i = np.minimum(np.floor(u).astype(int), dep.size - 2)
    f = u - i
    return disp[b, i] * (1 - f) + disp[b, i + 1] * f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--probe", type=int, required=True)
    ap.add_argument("--shank", type=int, default=0)
    ap.add_argument("--t0", type=float, required=True, help="stretch to use, viewer s")
    ap.add_argument("--t1", type=float, required=True)
    ap.add_argument("--data", default=str(Path(__file__).resolve().parent.parent / "data"))
    ap.add_argument("--alf", default="D:/temp/fd016_share/alf")
    a = ap.parse_args()
    from scipy import signal

    raw = Path(a.raw)
    man = json.load(open(raw / "manifest.json"))
    pr = man["probes"][a.probe]
    sh = pr["shanks"][a.shank]
    fs, tile, v = man["fs"], man["tile"], man["variant"]
    k = np.fromfile(raw / pr["clock"], "<f8").reshape(-1, 2)
    cs, ct = k[:, 0], k[:, 1]
    lin = np.polyfit(ct, cs, 1)                       # straight-line clock, for comparison
    disp = np.fromfile(raw / pr["motion"]["file"], "<f4").reshape(pr["motion"]["shape"])[a.shank]
    T0, T1 = int(np.interp(a.t0, ct, cs) // tile), int(np.interp(a.t1, ct, cs) // tile)
    s0 = T0 * tile
    x = np.concatenate([np.cumsum(np.frombuffer(zlib.decompress(
        (raw / f"p{a.probe}/sh{a.shank}/0.{T}").read_bytes()), "<i2"), dtype=np.int16).reshape(384, tile)
        for T in range(T0, T1 + 1)], axis=1)
    n = x.shape[1]
    sos = np.asarray(man["highpass"]["sos"])
    y_site = man["row_pitch_um"] * (np.arange(384) // 2)
    uv = sh["uv_per_count"]

    meta = json.load(open(Path(a.data) / "meta.json"))
    rows = meta["rows"]
    st = np.load(Path(a.alf) / "spikes.times.npy")
    sc = np.load(Path(a.alf) / "spikes.clusters.npy")
    t0, t1 = np.interp([s0 + HALF + 1, s0 + n - HALF - 2], cs, ct)
    i0, i1 = np.searchsorted(st, [t0, t1])
    st, sc = st[i0:i1], sc[i0:i1]

    units = [i for i in range(meta["n_rows"]) if rows["probe"][i] == a.probe and rows["shank"][i] == a.shank]
    spikes = {}
    for i in units:
        t = st[sc == rows["cluster"][i]]    # the viewer's rows are in depth order; "cluster" is the ALF index
        if t.size >= 30:
            spikes[i] = t
    print(f"probe {a.probe} shank {a.shank}: {len(units)} QC units, {len(spikes)} with 30 or more spikes in the stretch")

    sta = {i: {} for i in spikes}
    sta_lin = {i: {} for i in spikes}
    idx = np.arange(-HALF, HALF + 1)
    for c in range(384):
        near = [i for i in spikes if abs(rows["y_um"][i] - y_site[c]) <= 150]
        if not near or c in sh["dead"]:
            continue
        vv = signal.sosfiltfilt(sos, x[c].astype(np.float64) * v["scale"] + v["offset"]) * uv
        for i in near:
            s = np.rint(np.interp(spikes[i], ct, cs) - s0).astype(int)
            sta[i][c] = vv[s[:, None] + idx[None, :]].mean(0)
            s2 = np.clip(np.rint(np.polyval(lin, spikes[i]) - s0).astype(int), HALF, n - HALF - 1)
            sta_lin[i][c] = vv[s2[:, None] + idx[None, :]].mean(0)

    out = []
    for i in spikes:
        chans = sorted(sta[i])
        W = np.array([sta[i][c] for c in chans])               # sites x samples
        W2 = np.array([sta_lin[i][c] for c in chans])
        amp = -W.min(axis=1)
        kk = int(np.argmax(amp))
        lag = int(np.argmin(W[kk])) - HALF
        w = np.maximum(amp, 0) ** 2
        yc = float((w * y_site[chans]).sum() / w.sum())
        d = float(displacement(pr["motion"], disp, np.interp(spikes[i], ct, cs), rows["y_um"][i], fs).mean())
        out.append(dict(row=i, cluster=rows["cluster"][i], n_spikes=int(spikes[i].size),
                        amp_uv=float(amp[kk]), lag_samples=lag,
                        amp_straight_line_uv=float(-W2[kk].min()),
                        y_viewer_um=float(rows["y_um"][i]), y_raw_um=yc, motion_um=d,
                        y_peak_site_um=float(y_site[chans[kk]])))
    big = [o for o in out if o["amp_uv"] >= 30]
    lag = np.array([o["lag_samples"] for o in big])
    ratio = np.array([o["amp_straight_line_uv"] / o["amp_uv"] for o in big])
    d = np.array([o["motion_um"] for o in big])
    e0 = np.array([o["y_raw_um"] - o["y_viewer_um"] for o in big])
    stat = lambda e: dict(median=float(np.median(e)), median_abs=float(np.median(np.abs(e))),
                          p90_abs=float(np.percentile(np.abs(e), 90)))
    moved = np.abs(d) >= 5
    summ = dict(probe=a.probe, shank=a.shank, stretch_s=[a.t0, a.t1], n_units=len(out), n_units_30uv=len(big),
                lag_samples=dict(median=float(np.median(lag)), p5=float(np.percentile(lag, 5)),
                                 p95=float(np.percentile(lag, 95)), within_2=float(np.mean(np.abs(lag) <= 2))),
                straight_line_amp_ratio=dict(median=float(np.median(ratio)), p5=float(np.percentile(ratio, 5))),
                motion_abs_um=dict(median=float(np.median(np.abs(d))), max=float(np.abs(d).max())),
                depth_error_um={"y_raw - d (viewer)": stat(e0 - d), "y_raw": stat(e0), "y_raw + d": stat(e0 + d)},
                n_units_moved_5um=int(moved.sum()))
    if moved.sum() >= 3:
        summ["depth_error_um_moved_5um"] = {"y_raw - d (viewer)": stat((e0 - d)[moved]), "y_raw": stat(e0[moved]),
                                             "y_raw + d": stat((e0 + d)[moved])}
        # slope of (y_raw - y_unit) against d over all units: +1 if the rule is right, 0 if no shift is needed
        summ["slope_of_offset_on_motion"] = float(np.polyfit(d, e0, 1)[0])
    json.dump(dict(summary=summ, units=out), open(raw / f"validate_raw_p{a.probe}_sh{a.shank}.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
