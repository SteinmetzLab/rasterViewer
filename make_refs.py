"""Write reference cell values for the viewer's self-test (?test=1).

For a set of views (fixed cases plus random ones), computes reference.Data.cell_values and
writes refs/views.json (next to this file) plus one float32 file per view (H x W, row 0 at
the top). Cases deliberately include fractional row bounds, cells straddling shank
boundaries, very deep zoom, very wide grids, and both color modes, in each y-axis layout
(rows = sorted by depth, depth = plotted by depth, rastermap). Views without a "layout" field
are rows views; the first 22 are the original (rows-only) set, unchanged.
"""
import json
from pathlib import Path
import numpy as np
from reference import Data, DEPTH_H
sys_path_bursts = Path(__file__).resolve().parent / "bursts"
import sys; sys.path.insert(0, str(sys_path_bursts))
from detect_bursts import by_unit, detect  # noqa: E402

BURST_PARAMS = (100.0, 4.0, 2)   # the viewer's default burst definition

OUT = Path(__file__).resolve().parent / "refs"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    D = Data()
    fs, TS, TE, N = D.fs, D.meta["t_start"], D.meta["t_end"], D.n
    seg_starts = [s["start"] for s in D.meta["segments"]]
    views = []

    def add(name, t0, t1, n0, n1, W, H, norm=False, eps=0.1, layout=None, spikes=None):
        v = dict(name=name, T0=t0 * fs, dT=(t1 - t0) * fs / W, n0=n0, n1=n1,
                 W=W, H=H, norm=norm, eps=eps)
        if spikes:
            v["spikes"] = spikes; v["burst_params"] = list(BURST_PARAMS)
        if layout:
            v["layout"] = layout
            if layout == "depth":
                v["h"] = DEPTH_H
        views.append(v)

    add("full_abs_1400x800", TS, TE, 0, N, 1400, 800)
    add("full_norm_1400x800", TS, TE, 0, N, 1400, 800, True, 0.1)
    add("full_abs_wide_3000x1500", TS, TE, 0, N, 3000, 1500)
    add("full_abs_block3_466x266", TS, TE, 0, N, 466, 266)
    add("60s_all_abs", 1102, 1162, 0, N, 1400, 800)
    add("4s_all_norm_eps0.5", 1102, 1106, 0, N, 1400, 800, True, 0.5)
    b = seg_starts[9]
    add("straddle_shank_boundary", 1500, 1530, b - 7.3, b + 6.9, 900, 9)
    add("deep_50ms_20rows", 1102, 1102.05, 1242.4, 1262.4, 1400, 800)
    hi = int(np.argmax(D.r0))   # highest-rate neuron, so the deep view is not empty
    add("deep_30ms_2rows_busiest", 1102.001, 1102.031, hi - 0.75, hi + 1.25, 1200, 700)
    add("late_window_large_ticks", TE - 2, TE, 3000.5, 3400.25, 1400, 800, True, 0.0)
    add("dT_1.5_ticks_70ms", 1102.0, 1102.0 + 1400 * 1.5 / fs, 1240, 1340, 1400, 300)
    add("dT_21_ticks_1s", 1200.0, 1201.0, 0, N, 1400, 800, True, 0.1)
    rng = np.random.default_rng(7)
    for i in range(10):
        span = 10 ** rng.uniform(-1.0, np.log10(TE - TS))
        t0 = rng.uniform(TS, TE - span)
        nspan = 10 ** rng.uniform(np.log10(3), np.log10(N))
        n0 = rng.uniform(0, N - nspan)
        W, H = int(rng.integers(200, 2500)), int(rng.integers(100, 1300))
        norm = bool(rng.integers(0, 2))
        add(f"random{i}", t0, t0 + span, n0, n0 + nspan, W, H, norm, float(rng.choice([0, 0.1, 1.0])))

    # plotted by depth (Y in um, to scale)
    seg_lo, seg_hi = D.layout("depth")[4:6]
    YD = D.extent("depth")
    add("depth_full_abs", TS, TE, 0, YD, 1400, 800, layout="depth")
    add("depth_full_norm", TS, TE, 0, YD, 1244, 651, True, 0.1, layout="depth")
    add("depth_one_shank_60s", 1102, 1162, seg_lo[10], seg_hi[10], 1400, 800, layout="depth")
    add("depth_straddle_shank", 1500, 1530, seg_lo[9] - 37.3, seg_lo[9] + 41.9, 900, 23,
        layout="depth")
    # deep zooms anchored on the busiest unit's band start, so cells hold spikes and partial
    # overlaps (hi, the busiest row, is defined above)
    lo_hi = D.layout("depth")[0][hi]
    add("depth_deep_30um_1s", 1102, 1103, lo_hi - 13.7, lo_hi + 16.3, 1200, 700,
        layout="depth")
    add("depth_2um_tiny_cells", 1102, 1104, lo_hi - 1.13, lo_hi + 0.87, 800, 600,
        True, 0.1, layout="depth")
    add("depth_whole_probe", 1000, 1300, seg_lo[4], seg_hi[7], 2000, 1200, layout="depth")
    # plotted in Rastermap order
    add("rastermap_full_abs", TS, TE, 0, N, 1400, 800, layout="rastermap")
    add("rastermap_full_norm", TS, TE, 0, N, 2488, 1190, True, 0.1, layout="rastermap")
    add("rastermap_200_neurons_10s", 1102, 1112, 2600.4, 2800.9, 1400, 800, layout="rastermap")
    add("rastermap_2_neurons_50ms", 1102, 1102.05, 1000.25, 1002.25, 1400, 700, layout="rastermap")
    for i in range(6):
        kind = ["depth", "rastermap"][i % 2]
        Y = D.extent(kind)
        span = 10 ** rng.uniform(-0.5, np.log10(TE - TS))
        t0 = rng.uniform(TS, TE - span)
        ysp = Y * 10 ** rng.uniform(-3.0, 0.0)
        y0 = rng.uniform(0, Y - ysp)
        W, H = int(rng.integers(200, 2500)), int(rng.integers(100, 1300))
        add(f"random_{kind}{i}", t0, t0 + span, y0, y0 + ysp, W, H, bool(rng.integers(0, 2)),
            float(rng.choice([0, 0.1, 1.0])), layout=kind)

    # burst / tonic views (default burst definition), and the detection checksums the viewer
    # compares its own detection with
    order, _ = by_unit(D.ticks, D.rows.astype(np.int64), N)
    b_u, o_u, _ = detect(D.ticks[order], D.rows[order], *BURST_PARAMS, fs=fs)
    is_burst = np.zeros(D.ticks.size, bool); is_burst[order] = b_u
    is_onset = np.zeros(D.ticks.size, bool); is_onset[order] = o_u
    masks = {"burst": is_burst, "tonic": ~is_burst, "onset": is_onset}
    json.dump(dict(params=list(BURST_PARAMS), n_burst=int(is_burst.sum()), n_onset=int(is_onset.sum()),
                   sum_index_burst=int(np.flatnonzero(is_burst).sum()),
                   sum_index_onset=int(np.flatnonzero(is_onset).sum())),
              open(OUT / "bursts.json", "w"), indent=1)
    p6 = [s for s in D.meta["segments"] if s["probe"] == 5 and s["shank"] == 0][0]
    add("burst_full_abs", TS, TE, 0, N, 1400, 800, spikes="burst")
    add("tonic_full_norm", TS, TE, 0, N, 1244, 651, True, 0.1, spikes="tonic")
    add("onset_P6S0_60s", 1050, 1110, p6["start"], p6["end"], 1400, 800, spikes="onset")
    add("overlay_P6S0_screenshot", 1110.13, 1111.3, p6["start"], p6["end"], 1512, 590,
        spikes="overlay")
    add("burst_depth_full_norm", TS, TE, 0, D.extent("depth"), 900, 700, True, 0.1,
        layout="depth", spikes="burst")
    add("tonic_rastermap_10s", 1046, 1056, 0, N, 1400, 800, layout="rastermap", spikes="tonic")

    for i, v in enumerate(views):
        def ref(spk):
            return D.cell_values(v["T0"], v["dT"], v["n0"], v["n1"], v["W"], v["H"], v["norm"],
                                 v["eps"], layout=v.get("layout", "rows"), h=v.get("h", DEPTH_H),
                                 mask=None if spk is None else masks[spk])
        spk = v.get("spikes")
        a = ref("burst" if spk == "overlay" else spk)
        v["file"] = f"view{i:02d}.f32"
        a.astype("<f4").tofile(OUT / v["file"])
        if spk == "overlay":            # second output: tonic
            v["file2"] = f"view{i:02d}_b.f32"
            ref("tonic").astype("<f4").tofile(OUT / v["file2"])
        print(f"{v['name']:28s} {v['W']}x{v['H']}  max {a.max():.4g}  mean {a.mean():.4g}  "
              f"no-neuron cells {int((a == -2).sum())}  no-sample cols {int((a[0] == -1).sum())}")
    json.dump(views, open(OUT / "views.json", "w"), indent=1)


if __name__ == "__main__":
    main()
