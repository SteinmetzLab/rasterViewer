"""Follow-up to alternatives.py: hybrids of a hierarchical-clustering tree (local structure)
with a global coordinate (large-scale layout). Same held-out scoring; appends to alternatives.json.

  hc_average, flipped by X : keep the average-linkage tree, and at every node put the child whose
                             leaves have the lower mean X first (X = PC1 loading, spectral
                             coordinate, or Rastermap position). O(n) after the linkage.
  refine to sigma 3        : the smoothing refinement run down to a 3-neuron width.
"""
import json
import time

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

from alternatives import OUT, full_order, kernel, refine, score, spectral
from run_rastermap import binned_counts
from uniqueness import zrows


def order_tree_by(Z, coord):
    n = len(coord)
    mean = np.concatenate([coord.astype(float), np.zeros(n - 1)])
    cnt = np.concatenate([np.ones(n), np.zeros(n - 1)])
    for i, (a, b) in enumerate(Z[:, :2].astype(int)):
        cnt[n + i] = cnt[a] + cnt[b]
        mean[n + i] = (mean[a] * cnt[a] + mean[b] * cnt[b]) / cnt[n + i]
    out, stack = [], [2 * n - 2]
    while stack:
        node = stack.pop()
        if node < n:
            out.append(node)
        else:
            a, b = Z[node - n, :2].astype(int)
            first, second = (a, b) if mean[a] <= mean[b] else (b, a)
            stack += [second, first]
    return np.array(out)


def main():
    from rastermap import Rastermap
    X, _, _ = binned_counts(0.1)
    half = X.shape[1] // 2
    Xa = X[:, :half].copy()
    Zb = zrows(X[:, half:])
    Cb = (Zb @ Zb.T / Zb.shape[1]).astype(np.float32)
    del Zb, X
    n = Xa.shape[0]
    W = kernel(n)
    res = json.load(open(OUT))

    def record(name, order, secs, note=""):
        assert np.array_equal(np.sort(order), np.arange(n)), name
        r = score(Cb, W, order)
        r.update(seconds=round(secs, 1), note=note)
        res[name] = r
        json.dump(res, open(OUT, "w"), indent=1)
        print(name, json.dumps(r), flush=True)

    m0 = Rastermap(random_state=0, verbose=False)
    m0.fit(Xa)
    good = np.asarray(m0.igood).ravel()
    U = m0.Usv[good].astype(np.float64)
    Un = U / np.linalg.norm(U, axis=1, keepdims=True)
    t = time.time()
    Z = linkage(pdist(Un, "cosine"), "average")
    t_link = time.time() - t
    sp = spectral(Un)
    sp_rank = np.empty(len(sp)); sp_rank[sp] = np.arange(len(sp))
    for label, coord in [("PC1 loading", U[:, 0]), ("spectral coordinate", sp_rank),
                         ("rastermap position", m0.embedding[good, 0])]:
        t = time.time()
        o = order_tree_by(Z, coord)
        record(f"hc_average, flipped by {label}", full_order(o, good), t_link + time.time() - t,
               "linkage + flipping; the coordinate's own cost is extra")
    to_sub = np.cumsum(good) - 1
    rm_sub = to_sub[np.asarray(m0.isort)[: good.sum()]]
    t = time.time()
    record("rastermap seed 0 + refine (fine, to sigma 3)",
           full_order(refine(rm_sub, Un, sigmas=(25, 12, 6, 3)), good), time.time() - t)
    t = time.time()
    record("pc1 + refine (to sigma 3)",
           full_order(refine(np.argsort(U[:, 0]), Un, sigmas=(200, 100, 50, 25, 12, 6, 3)), good),
           time.time() - t)


if __name__ == "__main__":
    main()
