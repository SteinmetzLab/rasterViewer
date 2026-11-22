"""Quick benchmark of alternative ways to order the neurons, all scored on held-out data.

Every method is fitted on the FIRST half of the session (100 ms bins) and scored on the SECOND
half. Scores (all from the held-out correlation matrix of z-scored 100 ms counts):
  S        : neuron-level Rastermap-style objective, mean over pairs of corr x (-log distance),
             kernel normalized to mean 1 (see uniqueness.py).
  win_b    : mean pairwise correlation within non-overlapping blocks of b adjacent neurons,
             b = 3, 10, 30, 100. This is what matters when a viewer cell averages b rows.

Methods:
  rastermap seeds 0-2 (defaults); rastermap with n_splits = 1, 2; n_PCs = 64
  tiebreak   : rastermap seed 0, ties on the upsampled grid broken by parabolic interpolation
  refine     : iterated "smooth along the order, re-place each neuron at its best-matching
               position", coarse to fine (my own; a 1-D self-organizing-map style pass)
  spectral   : Fiedler vector of the normalized Laplacian of the (clipped) similarity
  pc1        : sort by loading on the first PC
  hc_*       : hierarchical clustering (average, ward) with and without optimal leaf ordering
Similarities for my own methods are cosine similarities of the neurons' 200 PC loadings
(the same denoised space Rastermap works in, population mean projected out).
Writes alternatives.json.
"""
import io
import json
import logging
import time
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
from scipy.linalg import eigh
from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import pdist

from run_rastermap import binned_counts
from uniqueness import kernel, zrows

HERE = Path(__file__).resolve().parent
OUT = HERE / "alternatives.json"


def win_corr(Cp, b):
    n = Cp.shape[0] // b
    v = [(Cp[i * b:(i + 1) * b, i * b:(i + 1) * b].sum() - np.trace(Cp[i * b:(i + 1) * b, i * b:(i + 1) * b]))
         / (b * (b - 1)) for i in range(n)]
    return float(np.mean(v))


def score(C, W, order):
    Cp = C[np.ix_(order, order)]
    n = len(order)
    r = {"S": round(float((Cp * W).sum() / (n * (n - 1) / 2)), 5)}
    for b in (3, 10, 30, 100):
        r[f"win_{b}"] = round(win_corr(Cp, b), 4)
    return r


def full_order(sub_order, good):
    """Order over all neurons: good ones in sub_order, the rest appended."""
    idx = np.flatnonzero(good)
    return np.concatenate([idx[sub_order], np.flatnonzero(~good)])


def refine(order, Un, sigmas=(200, 100, 50, 25, 12, 6), iters=3):
    order = np.asarray(order).copy()
    n = len(order)
    for s in sigmas:
        for _ in range(iters):
            T = gaussian_filter1d(Un[order], s, axis=0, mode="nearest")
            T /= np.linalg.norm(T, axis=1, keepdims=True) + 1e-12
            cc = Un @ T.T
            i = cc.argmax(1)
            ii = np.clip(i, 1, n - 2)
            a, b, c = cc[np.arange(n), ii - 1], cc[np.arange(n), ii], cc[np.arange(n), ii + 1]
            den = a - 2 * b + c
            off = np.where(den < 0, 0.5 * (a - c) / np.where(den == 0, 1, den), 0.0)
            order = np.argsort(i + np.clip(off, -0.5, 0.5), kind="stable")
    return order


def spectral(Un):
    Wm = np.clip(Un @ Un.T, 0, None)
    np.fill_diagonal(Wm, 0)
    d = Wm.sum(1)
    dm = 1 / np.sqrt(d)
    L = np.eye(len(d)) - (Wm * dm[:, None]) * dm[None, :]
    w, v = eigh(L, subset_by_index=[1, 1])
    return np.argsort(v[:, 0] * dm)


def main():
    from rastermap import Rastermap

    X, meta, _ = binned_counts(0.1)
    half = X.shape[1] // 2
    Xa = X[:, :half].copy()
    Zb = zrows(X[:, half:])
    Cb = (Zb @ Zb.T / Zb.shape[1]).astype(np.float32)
    del Zb, X
    n = Xa.shape[0]
    W = kernel(n)
    res = json.load(open(OUT)) if OUT.exists() else {}  # resume: finished fits are kept

    def record(name, order, secs, note=""):
        assert np.array_equal(np.sort(order), np.arange(n)), name
        r = score(Cb, W, order)
        r.update(seconds=round(secs, 1), note=note)
        res[name] = r
        json.dump(res, open(OUT, "w"), indent=1)
        print(name, json.dumps(r), flush=True)

    rng = np.random.default_rng(0)
    record("random", rng.permutation(n), 0)
    record("depth (viewer rows)", np.arange(n), 0)

    # rastermap baselines; seed 0 also with its stage timings from the log
    buf = io.StringIO()
    h = logging.StreamHandler(buf)
    logging.getLogger("rastermap").addHandler(h)
    logging.getLogger("rastermap").setLevel(logging.INFO)
    models = {}
    for seed in (0, 1, 2):
        if seed > 0 and f"rastermap seed {seed}" in res:
            continue
        t = time.time()
        m = Rastermap(random_state=seed, verbose=(seed == 0))
        m.fit(Xa)
        record(f"rastermap seed {seed}", np.asarray(m.isort), time.time() - t)
        models[seed] = m
    res["_rastermap_seed0_log"] = buf.getvalue().splitlines()
    m0 = models[0]
    good = np.asarray(m0.igood).ravel()
    U = m0.Usv[good].astype(np.float64)
    Un = U / np.linalg.norm(U, axis=1, keepdims=True)

    for label, kw in [("rastermap n_splits=1", dict(n_splits=1)),
                      ("rastermap n_splits=2", dict(n_splits=2)),
                      ("rastermap n_PCs=64", dict(n_PCs=64)),
                      ("rastermap n_clusters=50", dict(n_clusters=50)),
                      ("rastermap n_clusters=150", dict(n_clusters=150))]:
        if label in res:
            continue
        t = time.time()
        try:
            mm = Rastermap(random_state=0, verbose=False, **kw)
            mm.fit(Xa)
            record(label, np.asarray(mm.isort), time.time() - t)
        except Exception as e:  # report, don't hide
            res[label] = {"error": repr(e)}
            print(label, "ERROR", repr(e), flush=True)

    # tie-break on rastermap's grid
    t = time.time()
    Xrec = m0.U_upsampled
    Uz = (U - U.mean(1, keepdims=True)) / U.std(1, keepdims=True)
    cc = np.maximum(0, Xrec @ Uz.T)  # grid x neurons, as in rastermap's upsampling
    i = cc.argmax(0)
    ng = cc.shape[0]
    ii = np.clip(i, 1, ng - 2)
    cols = np.arange(cc.shape[1])
    a, b, c = cc[ii - 1, cols], cc[ii, cols], cc[ii + 1, cols]
    den = a - 2 * b + c
    off = np.where(den < 0, 0.5 * (a - c) / np.where(den == 0, 1, den), 0.0)
    sub = np.argsort(i + np.clip(off, -0.5, 0.5), kind="stable")
    # my grid assignment must reproduce rastermap's embedding; its own isort breaks ties with
    # numpy's default (unstable) argsort, so tie order there is arbitrary
    assert np.allclose(i / (ng / m0.n_clusters), m0.embedding[good, 0]), "grid mismatch"
    to_sub = np.cumsum(good) - 1
    rm_sub = to_sub[np.asarray(m0.isort)[: good.sum()]]
    assert good[np.asarray(m0.isort)[: good.sum()]].all()
    record("rastermap seed 0 + tiebreak", full_order(sub, good), time.time() - t,
           "extra seconds on top of the fit")

    t = time.time()
    record("rastermap seed 0 + refine", full_order(refine(rm_sub, Un), good), time.time() - t,
           "extra seconds on top of the fit")
    t = time.time()
    record("rastermap seed 0 + refine (fine only)",
           full_order(refine(rm_sub, Un, sigmas=(25, 12, 6)), good), time.time() - t,
           "extra seconds on top of the fit")

    t = time.time()
    pc1 = np.argsort(U[:, 0])
    record("pc1 sort", full_order(pc1, good), time.time() - t, "needs the PCs (about 3 s)")
    t = time.time()
    record("pc1 + refine", full_order(refine(pc1, Un), good), time.time() - t)

    t = time.time()
    sp = spectral(Un)
    record("spectral", full_order(sp, good), time.time() - t, "needs the PCs (about 3 s)")
    t = time.time()
    record("spectral + refine", full_order(refine(sp, Un), good), time.time() - t,
           "extra seconds on top of spectral")

    y = pdist(Un, "cosine")
    for method in ("average", "ward"):
        t = time.time()
        Z = linkage(pdist(Un, "euclidean") if method == "ward" else y, method)
        record(f"hc_{method}", full_order(leaves_list(Z), good), time.time() - t)
        t = time.time()
        Zo = optimal_leaf_ordering(Z, y)
        olo = leaves_list(Zo)
        record(f"hc_{method} + optimal leaf order", full_order(olo, good), time.time() - t,
               "extra seconds on top of the linkage")
        if method == "ward":
            t = time.time()
            record("hc_ward + optimal leaf order + refine (fine only)",
                   full_order(refine(olo, Un, sigmas=(25, 12, 6)), good), time.time() - t)


if __name__ == "__main__":
    main()
