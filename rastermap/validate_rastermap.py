"""Checks on the Rastermap orderings written by run_rastermap.py.

1. Stability: rank correlation and local-neighborhood overlap between the main ordering and
   variants (other k-means seed, 25 ms and 250 ms bins, fit on the first half only).
2. Does it put similar neurons next to each other, on data it did not see? The ordering fitted
   on the first half is scored on the second half: mean Pearson correlation (100 ms bins)
   between each neuron and its neighbors k = 1..5 positions away, compared with the
   probe/shank/y order used by the viewer and with random orders.
3. Composition: how often adjacent neurons in the order come from different probes, and
   whether position tracks firing rate.
4. Figure: whole session at 1 s bins, rows in depth order vs Rastermap order.

Writes validation.json and rastermap_check.png in this folder.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from run_rastermap import binned_counts

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(0)


def load(tag):
    z = np.load(HERE / f"rastermap_{tag}.npz")
    return {k: z[k] for k in z.files}


def rank_of(isort):
    r = np.empty(isort.size, int)
    r[isort] = np.arange(isort.size)
    return r


def neighbor_overlap(ra, rb, k=20):
    """Mean fraction of each neuron's k nearest (by rank) in order a that are within k in b."""
    ia, ib = np.argsort(ra), np.argsort(rb)
    out = []
    for i in range(ra.size):
        na = set(ia[max(0, ra[i] - k // 2): ra[i] + k // 2 + 1]) - {i}
        nb = set(ib[max(0, rb[i] - k // 2): rb[i] + k // 2 + 1]) - {i}
        out.append(len(na & nb) / max(1, len(na)))
    return float(np.mean(out))


def adjacent_corr(Z, order, lags=(1, 2, 3, 5, 10, 50)):
    """Mean correlation between rows order[p] and order[p+lag] of row-z-scored Z."""
    Zs = Z[order]
    return {lag: float(np.nanmean(np.sum(Zs[:-lag] * Zs[lag:], axis=1) / Z.shape[1]))
            for lag in lags}


def main():
    meta = json.load(open(HERE.parent / "data" / "meta.json"))
    probe = np.asarray(meta["rows"]["probe"])
    r0 = np.asarray(meta["rows"]["r0"])
    main_ = load("main")
    rm = rank_of(main_["isort"])
    out = {"n": int(rm.size)}

    # 1. stability
    stab = {}
    for tag in ["seed1", "bin25ms", "bin250ms", "firsthalf"]:
        v = load(tag)
        rv = rank_of(v["isort"])
        rho = spearmanr(rm, rv).statistic
        if rho < 0:  # Rastermap's direction is arbitrary
            rv = rv.size - 1 - rv
        stab[tag] = dict(spearman_abs=round(abs(float(rho)), 3),
                         neighbor_overlap_k20=round(neighbor_overlap(rm, rv, 20), 3),
                         neighbor_overlap_k100=round(neighbor_overlap(rm, rv, 100), 3))
    out["stability_vs_main"] = stab
    out["neighbor_overlap_chance_k20"] = round(20 / rm.size, 4)

    # 2. held-out neighbor similarity (second half, 100 ms bins)
    X, _, _ = binned_counts(0.1)
    half = X.shape[1] // 2
    Z = X[:, half:].astype(np.float64)
    Z -= Z.mean(1, keepdims=True)
    sd = Z.std(1, keepdims=True)
    Z = np.divide(Z, sd, out=np.full_like(Z, np.nan), where=sd > 0)
    fh = load("firsthalf")["isort"]
    depth = np.arange(rm.size)  # viewer rows are already probe, shank, y order
    res = {"rastermap_fit_on_first_half": adjacent_corr(Z, fh),
           "rastermap_fit_on_all (not held out)": adjacent_corr(Z, main_["isort"]),
           "depth_order (viewer rows)": adjacent_corr(Z, depth)}
    rnd = [adjacent_corr(Z, rng.permutation(rm.size)) for _ in range(5)]
    res["random_order (mean of 5)"] = {k: float(np.mean([r[k] for r in rnd])) for k in rnd[0]}
    out["second_half_mean_corr_by_lag"] = {k: {l: round(c, 4) for l, c in v.items()}
                                           for k, v in res.items()}

    # 3. composition
    o = main_["isort"]
    out["adjacent_pairs_crossing_probe"] = dict(
        rastermap=round(float(np.mean(probe[o][1:] != probe[o][:-1])), 3),
        depth_order=round(float(np.mean(probe[1:] != probe[:-1])), 4))
    out["spearman_position_vs_log_rate"] = round(float(spearmanr(rm, np.log(r0 + 1e-3)).statistic), 3)
    # probe makeup of 10 equal blocks along the order
    blocks = np.array_split(o, 10)
    out["probe_counts_in_10_blocks"] = [np.bincount(probe[b], minlength=probe.max() + 1).tolist()
                                       for b in blocks]
    json.dump(out, open(HERE / "validation.json", "w"), indent=1)
    print(json.dumps(out, indent=1))

    # 4. figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False,
                         "axes.spines.right": False})
    X1, _, _ = binned_counts(1.0)
    Zs = X1 - X1.mean(1, keepdims=True)
    Zs /= np.maximum(Zs.std(1, keepdims=True), 1e-9)
    fig, axs = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    for ax, ordr, title in [(axs[0], depth, "Rows in probe, shank, depth order (viewer)"),
                            (axs[1], o, "Rows in Rastermap order (100 ms bins, whole session)")]:
        im = ax.imshow(Zs[ordr], aspect="auto", cmap="gray_r", vmin=0, vmax=2.5,
                       interpolation="none", extent=[meta["t_start"], meta["t_start"] + X1.shape[1],
                                                     rm.size, 0])
        ax.set_title(title, loc="left")
        ax.set_ylabel("Neuron (position in order)")
        ax.tick_params(axis="y", pad=16)
        # probe identity strip on the left edge
        ax2 = ax.inset_axes([-0.018, 0, 0.012, 1])
        ax2.imshow(probe[ordr][:, None], aspect="auto", cmap="tab10", vmin=0, vmax=9,
                   interpolation="none")
        ax2.axis("off")
    axs[1].set_xlabel("Time (s)")
    cb = fig.colorbar(im, ax=axs, shrink=0.5, pad=0.01)
    cb.set_label("Firing rate (z-score, 1 s bins)")
    fig.savefig(HERE / "rastermap_check.png", dpi=110, bbox_inches="tight")


if __name__ == "__main__":
    main()
