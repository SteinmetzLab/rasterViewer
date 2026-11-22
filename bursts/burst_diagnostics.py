"""Diagnostic figures for choosing a burst definition (writes figs/ and summary.json here).

    python burst_diagnostics.py

Figures:
  1 ISI return maps: each spike's preceding vs following ISI, for groups of units
  2 per-unit fraction of spikes in bursts vs position along the shank, per probe and criterion
  3 example raster, Probe 6 shank 0, 1110.1-1111.3 s, burst vs tonic spikes per criterion
  4 burst structure: sizes, and ISI by position within the burst
  5 time course of tonic and burst firing, Probe 6 shank 0 lower half
  6 sensitivity: % of spikes in bursts over a grid of silence x max-ISI, deep vs superficial

Unit groups use position along the shank only (no histology): "deep" = y < 1000 um above the
deepest site, "superficial" = y > 2000 um. Probe 6 shank 0 is split at 1400 um, where its
dense deep population ends (see the viewer's plotted-by-depth view).
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pack_data import load_spikes  # noqa: E402
from detect_bursts import CRITERIA, FS, by_unit, detect  # noqa: E402

FIG = HERE / "figs"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
C_BURST, C_TONIC = "#d6249f", "#2ca02c"


def main():
    meta = json.load(open(HERE.parent / "data" / "meta.json"))
    R = meta["rows"]; n_units = meta["n_rows"]
    probe = np.array(R["probe"]); shank = np.array(R["shank"]); y = np.array(R["y_um"])
    ticks, rows = load_spikes(HERE.parent / "data")
    order, off = by_unit(ticks, rows.astype(np.int64), n_units)
    t = ticks[order].astype(np.int64); r = rows[order].astype(np.int64)
    same_prev = np.r_[False, r[1:] == r[:-1]]
    same_next = np.r_[r[1:] == r[:-1], False]
    pre_ms = np.where(same_prev, np.r_[0, np.diff(t)], -1) * 1000.0 / FS
    post_ms = np.where(same_next, np.r_[np.diff(t), 0], -1) * 1000.0 / FS
    unit_y = y[r]
    p6s0 = (probe == 5) & (shank == 0)
    groups = {
        "Probe 6 shank 0, y < 1400 um": p6s0 & (y < 1400),
        "Probe 6 shank 0, y > 1400 um": p6s0 & (y >= 1400),
        "All shanks, y < 1000 um (deep)": y < 1000,
        "All shanks, y > 2000 um (superficial)": y > 2000,
    }
    gmask = {k: v[r] for k, v in groups.items()}
    res = {name: detect(t, r, *par) for name, par in CRITERIA.items()}
    summary = {"n_spikes": int(t.size), "groups": {k: int(v.sum()) for k, v in groups.items()},
               "criteria": {}}

    # ---- 1: return maps
    fig, axs = plt.subplots(1, 4, figsize=(15, 4), constrained_layout=True)
    bins = np.linspace(np.log10(0.3), np.log10(3e4), 120)
    for ax, (g, m) in zip(axs, gmask.items()):
        ok = m & (pre_ms > 0) & (post_ms > 0)
        h, xe, ye, im = ax.hist2d(np.log10(pre_ms[ok]), np.log10(post_ms[ok]), bins=[bins, bins],
                                  norm=LogNorm(), cmap="magma")
        for (S, M, K), ls in zip(CRITERIA.values(), ["-", "--", ":"]):
            x0 = np.log10(max(S, 0.3))
            ax.plot([x0, x0, bins[-1]], [bins[0], np.log10(M), np.log10(M)], color="c", ls=ls, lw=1)
        ax.set_title(f"{g}\n{groups[g].sum()} units, {ok.sum():,} spikes", fontsize=9)
        tk = [1, 10, 100, 1000, 10000]
        ax.set_xticks(np.log10(tk), ["1", "10", "100", "1k", "10k"])
        ax.set_yticks(np.log10(tk), ["1", "10", "100", "1k", "10k"])
        ax.set_xlabel("ISI before the spike (ms)"); ax.set_ylabel("ISI after the spike (ms)")
    fig.suptitle("ISI return maps. LTS bursts' first spikes sit at long 'before' and short 'after' "
                 "(lower right). Cyan: solid LTS 100 ms / 4 ms, dashed 50 ms / 6 ms, dotted no silence / 8 ms")
    fig.savefig(FIG / "1_isi_return_maps.png", dpi=130); plt.close(fig)

    # ---- per-unit burst fraction
    nsp = np.bincount(r, minlength=n_units)
    frac = {name: np.bincount(r, weights=b.astype(float), minlength=n_units) / np.maximum(nsp, 1)
            for name, (b, o, _) in res.items()}
    onsets = {name: np.bincount(r, weights=o.astype(float), minlength=n_units)
              for name, (b, o, _) in res.items()}

    # ---- 2: burst fraction vs depth
    fig, axs = plt.subplots(len(CRITERIA), 7, figsize=(16, 9), sharey=True, sharex="row",
                            constrained_layout=True)
    for i, name in enumerate(CRITERIA):
        for p in range(7):
            ax = axs[i, p]
            for s in range(4):
                m = (probe == p) & (shank == s)
                ax.scatter(100 * frac[name][m], y[m], s=4, alpha=0.7, label=f"S{s}")
            if i == 0: ax.set_title(f"Probe {p + 1}")
            if p == 0: ax.set_ylabel(f"{name.split(':')[0]}\ny (um above deepest site)")
            if i == len(CRITERIA) - 1: ax.set_xlabel("% of spikes in bursts")
        axs[i, 0].legend(fontsize=7, markerscale=2, frameon=False)
    fig.suptitle("Per-unit fraction of spikes in bursts vs position along the shank (shallow at top)")
    fig.savefig(FIG / "2_burst_fraction_vs_depth.png", dpi=110); plt.close(fig)

    # ---- 3: example raster, P6S0
    t0s, t1s = 1110.1, 1111.3
    units = np.flatnonzero(p6s0)
    fig, axs = plt.subplots(len(CRITERIA), 1, figsize=(14, 4.2 * len(CRITERIA)),
                            constrained_layout=True, sharex=True)
    for ax, (name, (b, o, _)) in zip(axs, res.items()):
        for u in units:
            a, z = off[u], off[u + 1]
            tt = t[a:z] / FS; sel = (tt >= t0s) & (tt < t1s)
            if not sel.any(): continue
            bb = b[a:z][sel]; oo = o[a:z][sel]; xx = tt[sel]
            ax.scatter(xx[~bb], np.full((~bb).sum(), y[u]), marker="|", s=18, c=C_TONIC, lw=0.8)
            ax.scatter(xx[bb], np.full(bb.sum(), y[u]), marker="|", s=18, c=C_BURST, lw=0.8)
            ax.scatter(xx[oo], np.full(oo.sum(), y[u]), marker="v", s=10, c="k", lw=0)
        nb = int(sum(o[off[u]:off[u + 1]][(t[off[u]:off[u + 1]] / FS >= t0s) &
                                          (t[off[u]:off[u + 1]] / FS < t1s)].sum() for u in units))
        ax.set_title(f"{name}: {nb} bursts in this window (magenta burst spikes, green tonic, "
                     "triangle = burst onset)", fontsize=9)
        ax.set_ylabel("y (um above deepest site)")
    axs[-1].set_xlabel("Time (s)")
    fig.suptitle("Probe 6 shank 0, all units at their position along the shank, 1110.1-1111.3 s "
                 "(the window of the viewer screenshot)")
    fig.savefig(FIG / "3_example_raster_P6S0.png", dpi=120); plt.close(fig)

    # ---- 4: burst structure
    fig, axs = plt.subplots(2, len(CRITERIA), figsize=(14, 7), constrained_layout=True)
    for i, (name, (b, o, bid)) in enumerate(res.items()):
        for g, col in [("Probe 6 shank 0, y < 1400 um", "k"), ("All shanks, y > 2000 um (superficial)", "tab:orange")]:
            m = gmask[g] & b
            ids, size = np.unique(bid[m], return_counts=True)
            axs[0, i].hist(size, bins=np.arange(1.5, 12.5), histtype="step", color=col,
                           density=True, label=f"{g} ({ids.size:,} bursts)")
            # ISI k within the burst (between spike k and k+1), for bursts of >= 3 spikes
            first = np.flatnonzero(o & gmask[g])
            big = first[size[np.searchsorted(ids, bid[first])] >= 3] if ids.size else first[:0]
            med, lo_, hi_ = [], [], []
            for k in range(1, 4):
                ok = big[(big + k) < t.size]
                ok = ok[bid[ok + k] == bid[ok]] if ok.size else ok
                v = (t[ok + k] - t[ok + k - 1]) * 1000.0 / FS
                med.append(np.median(v) if v.size else np.nan)
                lo_.append(np.percentile(v, 25) if v.size else np.nan)
                hi_.append(np.percentile(v, 75) if v.size else np.nan)
            axs[1, i].errorbar(np.arange(1, 4), med, yerr=[np.subtract(med, lo_), np.subtract(hi_, med)],
                               color=col, marker="o", capsize=3, label=g)
        axs[0, i].set_title(name.split(":")[0], fontsize=9); axs[0, i].set_xlabel("spikes per burst")
        axs[0, i].legend(fontsize=7, frameon=False)
        axs[1, i].set_xticks([1, 2, 3], ["1st", "2nd", "3rd"]); axs[1, i].set_xlabel("ISI within burst (bursts of >= 3 spikes)")
        axs[1, i].set_ylabel("ms (median, IQR)")
    fig.suptitle("Burst structure. LTS bursts typically lengthen ISI by ISI (decelerate); fast tonic trains do not")
    fig.savefig(FIG / "4_burst_structure.png", dpi=120); plt.close(fig)

    # ---- 5: time course, P6S0 deep
    TS, TE = meta["t_start"], meta["t_end"]
    edges = np.arange(np.floor(TS), np.ceil(TE) + 1)
    m = gmask["Probe 6 shank 0, y < 1400 um"]
    nU = groups["Probe 6 shank 0, y < 1400 um"].sum()
    name = list(CRITERIA)[0]
    b, o, _ = res[name]
    tonic = np.histogram(t[m & ~b] / FS, edges)[0] / nU
    bon = np.histogram(t[m & o] / FS, edges)[0] / nU
    fig, axs = plt.subplots(3, 1, figsize=(14, 8), constrained_layout=True)
    for ax, (a, z) in zip(axs, [(TS, TE), (1050, 1170), (1100, 1120)]):
        sel = (edges[:-1] >= a) & (edges[:-1] < z)
        if z - a < 50:   # finer bins for the short window
            e2 = np.arange(a, z + 0.1, 0.1)
            ax.plot(e2[:-1], np.histogram(t[m & ~b] / FS, e2)[0] / nU / 0.1, color=C_TONIC, lw=0.8, label="tonic spikes/s per unit")
            ax.plot(e2[:-1], np.histogram(t[m & o] / FS, e2)[0] / nU / 0.1 * 10, color=C_BURST, lw=0.8, label="burst onsets/s per unit x 10")
            ax.set_title(f"{a}-{z} s, 100 ms bins")
        else:
            ax.plot(edges[:-1][sel], tonic[sel], color=C_TONIC, lw=0.6, label="tonic spikes/s per unit")
            ax.plot(edges[:-1][sel], bon[sel] * 10, color=C_BURST, lw=0.6, label="burst onsets/s per unit x 10")
            ax.set_title(f"{a:.0f}-{z:.0f} s, 1 s bins")
        ax.legend(fontsize=7, frameon=False, loc="upper right"); ax.set_xlabel("Time (s)")
    fig.suptitle(f"Probe 6 shank 0, y < 1400 um ({nU} units): tonic vs burst firing over time, {name.split(':')[0]}")
    fig.savefig(FIG / "5_timecourse_P6S0_deep.png", dpi=120); plt.close(fig)

    # ---- 6: sensitivity grid
    Ss, Ms = [0, 25, 50, 100, 200], [3, 4, 6, 8]
    grid = {g: np.zeros((len(Ss), len(Ms))) for g in ["All shanks, y < 1000 um (deep)", "All shanks, y > 2000 um (superficial)", "Probe 6 shank 0, y < 1400 um"]}
    for i, S in enumerate(Ss):
        for j, M in enumerate(Ms):
            b, _, _ = detect(t, r, S, M, 2)
            for g in grid:
                grid[g][i, j] = 100 * b[gmask[g]].mean()
    fig, axs = plt.subplots(1, 4, figsize=(16, 3.8), constrained_layout=True)
    for ax, (g, v) in zip(axs, grid.items()):
        im = ax.imshow(v, cmap="viridis", aspect="auto")
        for i in range(len(Ss)):
            for j in range(len(Ms)):
                ax.text(j, i, f"{v[i, j]:.1f}", ha="center", va="center", color="w", fontsize=8)
        ax.set_title(g + "\n% of spikes in bursts (>= 2 spikes)", fontsize=9)
        ax.set_xticks(range(len(Ms)), [f"<= {M}" for M in Ms]); ax.set_xlabel("max ISI in burst (ms)")
        ax.set_yticks(range(len(Ss)), [f">= {S}" for S in Ss]); ax.set_ylabel("silence before (ms)")
    ratio = grid["All shanks, y < 1000 um (deep)"] / np.maximum(grid["All shanks, y > 2000 um (superficial)"], 1e-9)
    im = axs[3].imshow(ratio, cmap="magma", aspect="auto")
    for i in range(len(Ss)):
        for j in range(len(Ms)):
            axs[3].text(j, i, f"{ratio[i, j]:.1f}", ha="center", va="center", color="w", fontsize=8)
    axs[3].set_title("deep / superficial ratio\n(higher = more specific to deep units)", fontsize=9)
    axs[3].set_xticks(range(len(Ms)), [f"<= {M}" for M in Ms]); axs[3].set_xlabel("max ISI in burst (ms)")
    axs[3].set_yticks(range(len(Ss)), [f">= {S}" for S in Ss])
    fig.savefig(FIG / "6_sensitivity.png", dpi=120); plt.close(fig)

    for name, (b, o, _) in res.items():
        summary["criteria"][name] = {
            "pct_spikes_in_bursts": {g: round(100 * float(b[m].mean()), 2) for g, m in gmask.items()},
            "bursts_per_unit_per_min": {g: round(float(o[m].sum()) / groups[g].sum() / ((TE - TS) / 60), 3)
                                        for g, m in gmask.items()},
            "n_bursts_total": int(o.sum()),
            "units_with_>=5pct_burst_spikes": int((frac[name] >= 0.05).sum()),
        }
    summary["sensitivity_pct_spikes_in_bursts"] = {g: v.round(2).tolist() for g, v in grid.items()}
    summary["sensitivity_axes"] = {"silence_ms": Ss, "max_isi_ms": Ms, "min_spikes": 2}
    json.dump(summary, open(HERE / "summary.json", "w"), indent=1)
    print(json.dumps(summary["criteria"], indent=1))


if __name__ == "__main__":
    main()
