"""Does a faster, harder search of the cluster order find a much better Rastermap score?

Uses the Steinmetz Lab fork of rastermap (D:/Dropbox/code/rastermap, fastsort.anneal_sort).
Data: first half of the session at 100 ms bins for fitting, second half held out.

Part 1 (fixed clusters, k-means seed 0): Rastermap's own greedy search vs annealing, at several
        budgets; then 1,000 annealing runs to see how often the best value is reached.
Part 2 (k-means seeds): for each seed, clusters -> annealed order -> neuron order, scored with
        the neuron-level score S on the fitted half and the held-out half. Greedy for the
        first 20 seeds as a paired comparison.
Writes fast_search.json.

    python fast_search.py part1
    python fast_search.py part2 --n_seeds 100
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from run_rastermap import binned_counts
from uniqueness import kernel, zrows

HERE = Path(__file__).resolve().parent
OUT = HERE / "fast_search.json"


def load_out():
    return json.load(open(OUT)) if OUT.exists() else {}


def canon(order):
    """order up to reversal, as a tuple"""
    o = tuple(int(x) for x in order)
    return min(o, o[::-1])


def part1():
    from rastermap import Rastermap, fastsort
    from rastermap.cluster import scaled_kmeans
    from rastermap.sort import traveling_salesman

    X, _, _ = binned_counts(0.1)
    Xa = X[:, : X.shape[1] // 2].copy()
    del X
    m = Rastermap(n_clusters=100, n_PCs=200, verbose=False).fit(Xa)   # for the PCs (and numba warm-up)
    good = np.asarray(m.igood).ravel()
    U_nodes = scaled_kmeans(m.Usv[good], n_clusters=100, random_state=0)[0]
    cc = (U_nodes @ U_nodes.T).astype(np.float32)
    n = cc.shape[0]
    res = {"n_clusters": int(n)}

    t = time.time()
    cc_g, inds_g, BBt = traveling_salesman(cc.copy(), locality=0.0)[:3]
    t_greedy = time.time() - t
    F_greedy = fastsort.objective(cc, BBt, inds_g)
    res["greedy"] = dict(seconds=round(t_greedy, 2), objective=round(F_greedy, 3),
                         objective_start=round(fastsort.objective(cc, BBt), 3))
    print("greedy", res["greedy"], flush=True)
    # is the greedy result a local optimum of the full neighborhood?
    p2, gain, nmv = fastsort.polish(cc.astype(np.float64), BBt.astype(np.float64),
                                    inds_g.astype(np.int64))
    res["greedy"]["polish_gain"] = round(gain, 3)
    res["greedy"]["polish_moves"] = nmv

    fastsort.anneal_sort(cc, n_restarts=2, n_steps=1000)          # compile
    res["budgets"] = []
    for n_steps in (50_000, 200_000, 1_000_000):
        for n_restarts in (1, 16):
            t = time.time()
            _, inds, _, info = fastsort.anneal_sort(cc, n_restarts=n_restarts, n_steps=n_steps,
                                                    random_state=1)
            r = dict(n_steps=n_steps, n_restarts=n_restarts, seconds=round(time.time() - t, 2),
                     best=round(info["objective"], 3),
                     runs_min=round(float(info["run_objectives"].min()), 3),
                     runs_median=round(float(np.median(info["run_objectives"])), 3))
            res["budgets"].append(r)
            print(r, flush=True)
    # the default budget plus the exhaustive polish
    t = time.time()
    _, inds, _, info = fastsort.anneal_sort(cc, n_restarts=16, n_steps=200_000, random_state=1,
                                            do_polish=True)
    res["default_plus_polish"] = dict(seconds=round(time.time() - t, 2),
                                      objective=round(info["objective"], 3),
                                      polish_moves=info["polish_moves"])
    print("default + polish", res["default_plus_polish"], flush=True)

    # 1,000 independent runs at the default budget, and 200 at 10x the budget
    for label, n_steps, R in (("runs_default_budget", 200_000, 1000),
                              ("runs_10x_budget", 2_000_000, 160)):
        t = time.time()
        _, inds, _, info = fastsort.anneal_sort(cc, n_restarts=R, n_steps=n_steps,
                                                random_state=2, return_all=True)
        F = info["run_objectives"]
        orders = info["run_orders"]
        # polish the 5 best runs: where do they end, and at the same order?
        top5 = np.argsort(F)[::-1][:5]
        pol = [fastsort.polish(cc.astype(np.float64), BBt.astype(np.float64), orders[i].copy())
               for i in top5]
        polF = [round(fastsort.objective(cc, BBt, p[0]), 3) for p in pol]
        info["objective"] = max(polF)
        inds = pol[int(np.argmax(polF))][0]
        res[label + "_top5_polished"] = dict(
            before=[round(float(F[i]), 3) for i in top5], after=polF,
            moves=[p[2] for p in pol], n_distinct_orders=len({canon(p[0]) for p in pol}))
        top = F.max()
        near = F > top - 1e-3
        distinct_top = len({canon(o) for o in orders[near]})
        res[label] = dict(
            n_runs=R, n_steps=n_steps, seconds=round(time.time() - t, 1),
            best=round(float(top), 3), polished_best=round(info["objective"], 3),
            percentiles={str(q): round(float(np.percentile(F, q)), 3)
                         for q in (0, 1, 10, 50, 90, 99, 100)},
            n_runs_at_best=int(near.sum()), n_distinct_orders_at_best=distinct_top,
            n_within_0p1pct_of_best=int((F > top * 0.999).sum()),
            n_beating_greedy=int((F > F_greedy + 1e-3).sum()),
            best_order=[int(x) for x in inds])
        print(label, {k: v for k, v in res[label].items() if k != "best_order"}, flush=True)
    out = load_out()
    out["part1"] = res
    json.dump(out, open(OUT, "w"), indent=1)


def part2(n_seeds, n_greedy):
    from rastermap import Rastermap
    from alternatives import score

    X, _, _ = binned_counts(0.1)
    half = X.shape[1] // 2
    Xa = X[:, :half].copy()
    Za, Zb = zrows(Xa), zrows(X[:, half:])
    Ca = (Za @ Za.T / Za.shape[1]).astype(np.float32)
    Cb = (Zb @ Zb.T / Zb.shape[1]).astype(np.float32)
    del Za, Zb, X
    n = Xa.shape[0]
    W = kernel(n)
    m0 = Rastermap(verbose=False).fit(Xa)
    Usv, Vsv = m0.Usv, m0.Vsv
    out = load_out()
    res = out.get("part2", {"anneal": {}, "greedy": {}})

    def run(seed, sorting):
        t = time.time()
        m = Rastermap(random_state=seed, sorting=sorting, n_restarts=16, verbose=False)
        m.fit(Usv=Usv, Vsv=Vsv)
        secs = time.time() - t
        o = np.asarray(m.isort)
        r = dict(train=score(Ca, W, o), test=score(Cb, W, o), seconds=round(secs, 1),
                 cluster_objective=round(float((m.cc.astype(np.float64) * m.BBt).sum()), 3))
        return r

    # the class fitted from saved PCs must reproduce the fit from data (seed 0, greedy)
    chk = Rastermap(random_state=0, verbose=False)
    chk.fit(Usv=Usv, Vsv=Vsv)
    assert np.array_equal(np.asarray(chk.isort), np.asarray(m0.isort)), "PC-input fit differs"

    for seed in range(n_seeds):
        for sorting, limit in (("anneal", n_seeds), ("greedy", n_greedy)):
            if seed >= limit or str(seed) in res[sorting]:
                continue
            res[sorting][str(seed)] = run(seed, sorting)
            r = res[sorting][str(seed)]
            print(seed, sorting, r["cluster_objective"], r["train"]["S"], r["test"]["S"],
                  r["seconds"], flush=True)
        out["part2"] = res
        json.dump(out, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["part1", "part2"])
    ap.add_argument("--n_seeds", type=int, default=100)
    ap.add_argument("--n_greedy", type=int, default=20)
    a = ap.parse_args()
    part1() if a.part == "part1" else part2(a.n_seeds, a.n_greedy)
