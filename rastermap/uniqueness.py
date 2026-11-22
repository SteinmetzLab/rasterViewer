"""Why do Rastermap orderings differ between runs, and do the different ones score differently?

Rastermap has three stages (rastermap 1.0 source: rastermap.py, cluster.py, sort.py, upsample.py):
  1. scaled k-means of the neurons' PC loadings into 100 clusters (kmeans++ init from the seed),
  2. ordering of the 100 clusters by a greedy local search (move or flip one segment at a time)
     that maximizes sum_ij cc_ij * BBt_ij, BBt = -log(eps + |x_i - x_j|) normalized; it stops
     when no single move improves the score by more than 1e-3,
  3. each neuron is placed at the best-matching point of a 1,000-point grid interpolated between
     the ordered cluster centers. Neurons on the same grid point tie, and argsort breaks the
     tie by input order (here: probe, shank, y).

This script:
  A. fits 10 seeds on the FIRST half (100 ms bins) and scores every ordering on the first half
     (train) and the second half (test) with a neuron-level version of Rastermap's objective,
     S(order) = mean over pairs i<j of C[o_i, o_j] * w(|i - j| / N), w = the same -log kernel
     normalized to mean 1 over pairs, C = correlation of the 100 ms counts. Also the mean
     correlation of neighbors at k = 1, 10, 50. Depth order and random orders for reference.
  B. for seed 0's clusters, reruns stage 2 from 30 random starting orders, to see whether the
     cluster ordering search itself lands in different local optima.
  C. counts ties from stage 3.
Writes uniqueness.json.
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from run_rastermap import binned_counts

HERE = Path(__file__).resolve().parent
N_SEEDS = 10


def zrows(X):
    Z = X.astype(np.float64)
    Z -= Z.mean(1, keepdims=True)
    sd = Z.std(1, keepdims=True)
    return np.divide(Z, sd, out=np.zeros_like(Z), where=sd > 0)


def kernel(n, eps=1e-3):
    d = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :]) / n
    W = -np.log(eps + d)
    iu = np.triu_indices(n, 1)
    W = np.triu(W, 1)
    return (W / W[iu].mean()).astype(np.float32)


def scores(C, W, order):
    Cp = C[np.ix_(order, order)]
    n = len(order)
    S = float((Cp * W).sum() / (n * (n - 1) / 2))
    adj = {k: float(np.mean(np.diagonal(Cp, k))) for k in (1, 10, 50)}
    return dict(S=round(S, 5), **{f"corr_k{k}": round(v, 4) for k, v in adj.items()})


def main():
    from rastermap import Rastermap
    from rastermap.sort import traveling_salesman

    X, meta, _ = binned_counts(0.1)
    half = X.shape[1] // 2
    Xa, Xb = X[:, :half].copy(), X[:, half:].copy()
    Za, Zb = zrows(Xa), zrows(Xb)
    Ca = (Za @ Za.T / Za.shape[1]).astype(np.float32)
    Cb = (Zb @ Zb.T / Zb.shape[1]).astype(np.float32)
    del Za, Zb
    n = X.shape[0]
    W = kernel(n)
    out = {"definition": __doc__, "per_seed": {}, "reference": {}}
    rng = np.random.default_rng(0)

    def both(order):
        return {"train_first_half": scores(Ca, W, order), "test_second_half": scores(Cb, W, order)}

    out["reference"]["depth_order"] = both(np.arange(n))
    out["reference"]["random_order"] = [both(rng.permutation(n)) for _ in range(3)]

    orders, models = [], []
    for seed in range(N_SEEDS):
        m = Rastermap(n_clusters=100, n_PCs=200, locality=0.0, time_lag_window=0,
                      random_state=seed, verbose=False)
        m.fit(Xa)
        o = np.asarray(m.isort)
        emb = m.embedding[:, 0]
        good = np.isfinite(emb)
        n_unique = int(np.unique(emb[good]).size)
        # Rastermap's own (cluster-level) objective, as maximized in stage 2
        S_clust = float((m.cc.astype(np.float32) * m.BBt.astype(np.float32)).sum())
        r = both(o)
        r.update(cluster_objective=round(S_clust, 2), n_distinct_positions=n_unique,
                 mean_neurons_per_position=round(good.sum() / n_unique, 2))
        out["per_seed"][seed] = r
        orders.append(o)
        models.append(m)
        print(seed, json.dumps(r), flush=True)

    # pairwise agreement between seeds
    ranks = []
    for o in orders:
        rk = np.empty(n, int)
        rk[o] = np.arange(n)
        ranks.append(rk)
    rho = np.array([[abs(spearmanr(a, b).statistic) for b in ranks] for a in ranks])
    out["seed_pairwise_abs_spearman"] = {
        "mean": round(float(rho[np.triu_indices(N_SEEDS, 1)].mean()), 3),
        "min": round(float(rho[np.triu_indices(N_SEEDS, 1)].min()), 3),
        "max": round(float(rho[np.triu_indices(N_SEEDS, 1)].max()), 3)}

    # B. stage 2 restarts on fixed clusters (seed 0)
    m = models[0]
    Un = m.U_nodes  # already in the sorted order
    cc0 = (Un @ Un.T).astype(np.float32)
    base = traveling_salesman(cc0.copy(), locality=0.0)
    S_base = float((base[0] * base[2]).sum())
    restarts = []
    for _ in range(30):
        p = rng.permutation(cc0.shape[0])
        res = traveling_salesman(cc0[np.ix_(p, p)].copy(), locality=0.0)
        inds = p[res[1]]
        S = float((res[0] * res[2]).sum())
        # agreement of the cluster order with the default one (direction-free)
        rr = abs(spearmanr(np.argsort(inds), np.arange(len(inds))).statistic)
        restarts.append(dict(S=round(S, 2), spearman_vs_default=round(float(rr), 3)))
    Ss = np.array([r["S"] for r in restarts])
    out["stage2_restarts_seed0"] = dict(
        default_start_S=round(S_base, 2), restart_S_min=float(Ss.min()),
        restart_S_median=float(np.median(Ss)), restart_S_max=float(Ss.max()),
        n_restarts_beating_default=int((Ss > S_base + 1e-3).sum()), restarts=restarts)
    json.dump(out, open(HERE / "uniqueness.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("definition", "per_seed")},
                     indent=1)[:4000])


if __name__ == "__main__":
    main()
