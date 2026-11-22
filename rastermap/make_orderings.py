"""Final orderings of the whole session for the viewer (uses the lab's rastermap fork).

    python make_orderings.py [--hp_period 100]

Fits Rastermap (defaults, seed 0, 100 ms bins) on the whole session and derives:
  rm_refine   Rastermap order with the local order refined (postsort.refine_isort)
  rm_tree     hierarchical-clustering tree, branches flipped to follow the Rastermap position
              (postsort.tree_isort)
and the same three after removing slow drift from the input (every cosine mode with period
longer than --hp_period seconds, see drift.py):
  hp_rastermap, hp_refine, hp_tree

Writes orderings.json here (all six, with isort / rank / parameters) and a check figure
orderings_check.png. The plain Rastermap fit must reproduce rastermap_main.json.
The spikes shown by the viewer are never filtered; only the sorting input is.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from drift import BIN, highpass
from run_rastermap import binned_counts

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hp_period", type=float, default=100.0)
    a = ap.parse_args()
    from rastermap import Rastermap, postsort

    X, meta, n_used = binned_counts(BIN)
    n = X.shape[0]
    out = {}

    def fit(Xin, tag, extra):
        m = Rastermap(n_clusters=100, n_PCs=200, locality=0.0, time_lag_window=0,
                      random_state=0, verbose=False).fit(Xin)
        orders = {"rastermap": np.asarray(m.isort), "refine": postsort.refine_isort(m),
                  "tree": postsort.tree_isort(m)}
        for name, o in orders.items():
            o = np.asarray(o).astype(int)
            assert np.array_equal(np.sort(o), np.arange(n)), (tag, name)
            rank = np.empty(n, int)
            rank[o] = np.arange(n)
            out[f"{tag}_{name}"] = dict(isort=o.tolist(), rank=rank.tolist(), **extra)
        return m

    base = dict(bin_ms=BIN * 1000, n_clusters=100, n_PCs=200, locality=0.0, time_lag_window=0,
                random_state=0, sorting="greedy", rastermap="SteinmetzLab fork of v1.0")
    fit(X, "rm", dict(params=dict(base, highpass_period_s=None)))
    main_isort = json.load(open(HERE / "rastermap_main.json"))["isort"]
    assert out["rm_rastermap"]["isort"] == main_isort, "plain fit differs from rastermap_main.json"
    Xh = highpass(X, a.hp_period)
    fit(Xh, "hp", dict(params=dict(base, highpass_period_s=a.hp_period)))

    notes = {
        "rm_rastermap": "Rastermap defaults (identical to rastermap_main.json)",
        "rm_refine": "Rastermap, then local refinement (smoothing widths 25, 12, 6, 3 neurons)",
        "rm_tree": "Average-linkage tree of the PC loadings, branches flipped to follow Rastermap",
    }
    for k in list(out):
        kind = k.split("_", 1)[1]
        out[k]["description"] = notes["rm_" + kind] + (
            f"; sorting input high-passed (cosine modes with period > {a.hp_period:g} s removed)"
            if k.startswith("hp_") else "")
    json.dump(dict(description="Orderings of the viewer's display rows. isort[k] is the display "
                   "row at position k; rank is the inverse.", cluster_id=meta["rows"]["cluster"],
                   orderings=out), open(HERE / "orderings.json", "w"))

    # check figure: whole session, 10 s bins, each neuron divided by its mean, as in the viewer's
    # normalized mode
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "Arial", "axes.spines.top": False,
                         "axes.spines.right": False})
    nb = X.shape[1] // 100
    R = X[:, : nb * 100].reshape(n, nb, 100).sum(2) / (100 * BIN)
    R = R / (R.mean(1, keepdims=True) + 0.1)
    keys = ["rm_rastermap", "rm_refine", "rm_tree", "hp_rastermap", "hp_refine", "hp_tree"]
    fig, axs = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
    for ax, k in zip(axs.ravel(), keys):
        o = np.array(out[k]["isort"])
        # average blocks of 5 neurons, as a viewer cell would
        img = R[o][: n // 5 * 5].reshape(n // 5, 5, nb).mean(1)
        im = ax.imshow(np.log10(np.maximum(img, 1e-3)), aspect="auto", cmap="viridis", vmin=-1,
                       vmax=1, interpolation="none",
                       extent=[meta["t_start"], meta["t_start"] + nb * 10, n, 0])
        kind = {"rastermap": "Rastermap", "refine": "Rastermap, locally refined",
                "tree": "Clustering tree, Rastermap layout"}[k.split("_", 1)[1]]
        ax.set_title(kind + (f"; sorted on periods under {a.hp_period:g} s only"
                             if k.startswith("hp_") else "; sorted on all timescales"),
                     loc="left", fontsize=10)
    for ax in axs[-1]:
        ax.set_xlabel("Time (s)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Neuron (position in order)")
    cb = fig.colorbar(im, ax=axs, shrink=0.5, pad=0.01)
    cb.set_label("Rate / (session mean + 0.1), log10")
    fig.savefig(HERE / "orderings_check.png", dpi=100, bbox_inches="tight")


if __name__ == "__main__":
    main()
