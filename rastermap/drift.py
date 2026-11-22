"""Slow drift, nonsense correlations, and what a high-pass filter does to the ordering.

Time courses are decomposed with the DCT-II (cosine modes k = 0, 1, 2, ...; mode k has period
2T/k, T = duration). Unlike the FFT it represents a trend across the session without leakage:
k = 1 is a monotone drift (down or up), k = 2 a peak or dip in the middle. Each mode is one real
number per neuron, i.e. one degree of freedom.

Part A  (whole session, 100 ms bins, z-scored, population mean projected out as Rastermap does)
        per band of periods: share of each neuron's variance, share of the variance in the
        200 PCs Rastermap sorts on, spread (SD) of pairwise correlations within the band, and
        the same spread under a null where every neuron keeps its own spectrum but the signs of
        its modes are flipped at random (independent neurons with the same slow statistics).
Part B  the same bands computed separately in the two halves: is a pair's within-band
        correlation in the first half reproduced in the second half? (reliability = Pearson r
        across the 9.3 M pairs)
Part C  Rastermap (and the local re-orderings) fitted on the first half after removing periods
        longer than a cutoff, scored on the second half.

    python drift.py A | B | C
Writes drift.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.fft import dct, idct

from run_rastermap import binned_counts
from uniqueness import kernel

HERE = Path(__file__).resolve().parent
OUT = HERE / "drift.json"
BIN = 0.1
K_EDGES = [1, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 4097]


def load_out():
    return json.load(open(OUT)) if OUT.exists() else {}


def save(key, val):
    out = load_out()
    out[key] = val
    json.dump(out, open(OUT, "w"), indent=1)


def zscore_rows(X):
    Z = X.astype(np.float32)
    Z -= Z.mean(1, keepdims=True)
    sd = Z.std(1, keepdims=True)
    return np.divide(Z, sd, out=np.zeros_like(Z), where=sd > 0)


def project_out_mean(Z):
    m = Z.mean(0, keepdims=True).T
    m /= np.sqrt((m**2).sum())
    return Z - (Z @ m) @ m.T


def highpass(X, period_s, bin_s=BIN):
    """Remove the mean and every cosine mode with period longer than period_s (per row)."""
    if period_s is None:
        return X.astype(np.float32)
    D = dct(X.astype(np.float32), type=2, norm="ortho", axis=1, workers=-1)
    T = X.shape[1] * bin_s
    kc = int(np.floor(2 * T / period_s)) + 1        # modes 0..kc-1 have period > period_s
    D[:, :kc] = 0
    return idct(D, type=2, norm="ortho", axis=1, workers=-1)


def bands(n_modes):
    edges = [e for e in K_EDGES if e < n_modes] + [n_modes]
    return list(zip(edges[:-1], edges[1:]))


def band_corr(A):
    """Correlation matrix of rows of A (mode coefficients in one band; no mean removal: the
    modes k >= 1 are already orthogonal to the mean)."""
    nrm = np.sqrt((A.astype(np.float64)**2).sum(1, keepdims=True))
    An = (A / np.where(nrm > 0, nrm, 1)).astype(np.float32)
    return An @ An.T


def offdiag(C):
    return C[np.triu_indices(C.shape[0], 1)]


def part_A():
    from rastermap import Rastermap
    X, _, _ = binned_counts(BIN)
    T = X.shape[1] * BIN
    Z = project_out_mean(zscore_rows(X))
    m = Rastermap(verbose=False).fit(X)
    del X
    good = np.asarray(m.igood).ravel()
    U = m.Usv[good] / m.sv                     # orthonormal neuron-space PCs (4325 x 200)
    D = dct(Z, type=2, norm="ortho", axis=1, workers=-1)
    tot = float((D[:, 1:].astype(np.float64)**2).sum())
    UD = U.T.astype(np.float32) @ D[good]      # PCs x modes
    tot_pc = float((UD[:, 1:].astype(np.float64)**2).sum())
    tot_pc10 = float((UD[:10, 1:].astype(np.float64)**2).sum())
    rng = np.random.default_rng(0)
    rows = []
    for k0, k1 in bands(D.shape[1]):
        A = D[:, k0:k1]
        C = band_corr(A)
        r = offdiag(C)
        As = A * rng.choice(np.float32([-1, 1]), size=A.shape)
        rs = offdiag(band_corr(As))
        ev = np.linalg.eigvalsh(C.astype(np.float64))[::-1]
        evs = np.linalg.eigvalsh(band_corr(As).astype(np.float64))[::-1]
        row = dict(k0=k0, k1=k1, dof=k1 - k0,
                   period_long_s=round(2 * T / k0, 1), period_short_s=round(2 * T / max(k1 - 1, 1), 1),
                   var_share=round(float((A.astype(np.float64)**2).sum()) / tot, 4),
                   var_share_200pc=round(float((UD[:, k0:k1].astype(np.float64)**2).sum()) / tot_pc, 4),
                   var_share_10pc=round(float((UD[:10, k0:k1].astype(np.float64)**2).sum()) / tot_pc10, 4),
                   sd_r=round(float(r.std()), 4), sd_r_null=round(float(rs.std()), 4),
                   mean_abs_r=round(float(np.abs(r).mean()), 4),
                   top_eig_share=round(float(ev[0] / ev.sum()), 4),
                   top_eig_share_null=round(float(evs[0] / evs.sum()), 4))
        rows.append(row)
        print(row, flush=True)
    save("A_bands_whole_session", rows)


def part_B():
    X, _, _ = binned_counts(BIN)
    half = X.shape[1] // 2
    T = half * BIN
    Ds = []
    for sl in (slice(0, half), slice(half, 2 * half)):
        Z = project_out_mean(zscore_rows(X[:, sl]))
        Ds.append(dct(Z, type=2, norm="ortho", axis=1, workers=-1))
    del X
    rng = np.random.default_rng(0)
    rows = []
    for k0, k1 in bands(half):
        C1, C2 = band_corr(Ds[0][:, k0:k1]), band_corr(Ds[1][:, k0:k1])
        r1, r2 = offdiag(C1), offdiag(C2)
        A = Ds[0][:, k0:k1]
        rs = offdiag(band_corr(A * rng.choice(np.float32([-1, 1]), size=A.shape)))
        rel = float(np.corrcoef(r1, r2)[0, 1])
        row = dict(k0=k0, k1=k1, dof_per_half=k1 - k0,
                   period_long_s=round(2 * T / k0, 1), period_short_s=round(2 * T / max(k1 - 1, 1), 1),
                   reliability=round(rel, 4), sd_r_half1=round(float(r1.std()), 4),
                   sd_r_half2=round(float(r2.std()), 4), sd_r_null=round(float(rs.std()), 4))
        rows.append(row)
        print(row, flush=True)
    save("B_split_half_reliability", rows)


def part_C(cutoffs=(None, 1000, 300, 100, 30, 10), evals=(None, 300, 100, 30)):
    from rastermap import Rastermap, postsort
    from alternatives import score

    X, _, _ = binned_counts(BIN)
    half = X.shape[1] // 2
    Xa, Xb = X[:, :half].copy(), X[:, half:2 * half].copy()
    del X
    n = Xa.shape[0]
    W = kernel(n)
    Ceval = {}
    for p in evals:                                   # held-out correlation matrices
        Zb = zscore_rows(highpass(Xb, p)).astype(np.float64)
        Ceval[str(p)] = (Zb @ Zb.T / Zb.shape[1]).astype(np.float32)
    # slow part only (periods > 300 s), of each half, to see whether an order groups drifts
    Cslow = {}
    for name, Xh in (("fit", Xa), ("test", Xb)):
        Zs = zscore_rows(Xh.astype(np.float32) - highpass(Xh, 300) - Xh.mean(1, keepdims=True))
        Cslow[name] = (Zs.astype(np.float64) @ Zs.T.astype(np.float64) / Zs.shape[1]).astype(np.float32)
    res = load_out().get("C_cutoffs", {})
    for p in cutoffs:
        t = time.time()
        Xf = highpass(Xa, p)
        m = Rastermap(verbose=False, sorting="anneal").fit(Xf)
        orders = {"rastermap": np.asarray(m.isort), "refine": postsort.refine_isort(m),
                  "tree": postsort.tree_isort(m)}
        for name, o in orders.items():
            r = {f"test_hp{e}": score(Ceval[str(e)], W, o) for e in evals}
            r["slow_fit_win30"] = score(Cslow["fit"], W, o)["win_30"]
            r["slow_test_win30"] = score(Cslow["test"], W, o)["win_30"]
            res[f"fit_hp{p}|{name}"] = r
            print(p, name, json.dumps(r), flush=True)
        print("cutoff", p, round(time.time() - t, 1), "s", flush=True)
        save("C_cutoffs", res)
    for name, o in (("depth", np.arange(n)), ("random", np.random.default_rng(0).permutation(n))):
        r = {f"test_hp{e}": score(Ceval[str(e)], W, o) for e in evals}
        r["slow_fit_win30"] = score(Cslow["fit"], W, o)["win_30"]
        r["slow_test_win30"] = score(Cslow["test"], W, o)["win_30"]
        res[name] = r
    save("C_cutoffs", res)


if __name__ == "__main__":
    {"A": part_A, "B": part_B, "C": part_C}[sys.argv[1]]()
