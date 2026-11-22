"""Rastermap ordering of the FD_016 2026-11-13/1 population (4,325 units).

Reads the viewer's spike files (../data: the packed spikes and meta.json), bins the
spikes into a neurons x time count matrix, and fits Rastermap (MouseLand, v1.0).

Usage (run with the rastermap venv, D:/temp/rastermap_env/Scripts/python.exe):
    python run_rastermap.py scale --bin_ms 100 --frac 0.1        # one timing run on a subset
    python run_rastermap.py full  --bin_ms 100 --tag main         # the full session

Each run prints one JSON line with wall time and peak memory (peak working set, Windows).
A `full` run writes rastermap_<tag>.npz and rastermap_<tag>.json in this folder.

Row indices are the viewer's display rows (0..4324, ordered probe, shank, y). Each row's
ALF cluster id is meta['rows']['cluster'][row].
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import psutil

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
import sys; sys.path.insert(0, str(HERE.parent))
from pack_data import load_spikes  # noqa: E402


def binned_counts(bin_s):
    """Spike counts per (row, bin) over [t_start, t_end), float32, shape (n_rows, n_bins)."""
    meta = json.load(open(DATA / "meta.json"))
    fs, n = meta["fs"], meta["n_rows"]
    t0, t1 = meta["t_start"], meta["t_end"]
    ticks, rows = load_spikes(DATA)
    assert ticks.size == rows.size == meta["n_spikes"]
    nb = int(np.floor((t1 - t0) / bin_s))
    b = np.floor((ticks.astype(np.float64) / fs - t0) / bin_s).astype(np.int64)
    ok = (b >= 0) & (b < nb)
    flat = rows[ok].astype(np.int64) * nb + b[ok]
    X = np.bincount(flat, minlength=n * nb).astype(np.float32).reshape(n, nb)
    return X, meta, int(ok.sum())


def peak_mb():
    return psutil.Process().memory_info().peak_wset / 2**20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["scale", "full"])
    ap.add_argument("--bin_ms", type=float, default=100.0)
    ap.add_argument("--frac", type=float, default=1.0, help="fraction of the session (from the start)")
    ap.add_argument("--n_clusters", type=int, default=100)
    ap.add_argument("--n_PCs", type=int, default=200)
    ap.add_argument("--locality", type=float, default=0.0)
    ap.add_argument("--time_lag_window", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0, help="Rastermap random_state (k-means init)")
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()

    from rastermap import Rastermap  # imported here so binning time/memory is separate

    tic = time.time()
    X, meta, n_used = binned_counts(a.bin_ms / 1000)
    if a.frac < 1:
        X = np.ascontiguousarray(X[:, : int(round(X.shape[1] * a.frac))])
    t_bin = time.time() - tic
    mem_bin = peak_mb()

    tic = time.time()
    model = Rastermap(n_clusters=a.n_clusters, n_PCs=a.n_PCs, locality=a.locality,
                      time_lag_window=a.time_lag_window, random_state=a.seed,
                      verbose=False)
    model.fit(X)
    t_fit = time.time() - tic

    rep = dict(mode=a.mode, bin_ms=a.bin_ms, frac=a.frac, n_rows=X.shape[0], n_bins=X.shape[1],
               X_mb=round(X.nbytes / 2**20), t_bin_s=round(t_bin, 1), t_fit_s=round(t_fit, 1),
               peak_mb_after_binning=round(mem_bin), peak_mb=round(peak_mb()),
               n_good=int(np.isfinite(model.embedding[:, 0]).sum()))
    print(json.dumps(rep), flush=True)

    if a.mode == "full":
        isort = np.asarray(model.isort)
        emb = np.asarray(model.embedding[:, 0])
        cl = np.asarray(meta["rows"]["cluster"])
        params = dict(n_clusters=a.n_clusters, n_PCs=a.n_PCs, locality=a.locality,
                      time_lag_window=a.time_lag_window, grid_upsample=10, normalize=True,
                      mean_time=True, bin_ms=a.bin_ms, t_start=meta["t_start"],
                      t_end_used=meta["t_start"] + X.shape[1] * a.bin_ms / 1000,
                      rastermap_version="1.0", random_state=a.seed, spikes_binned=n_used)
        np.savez(HERE / f"rastermap_{a.tag}.npz", isort=isort, embedding=emb,
                 embedding_clust=np.asarray(model.embedding_clust), cluster_id=cl,
                 cc=np.asarray(model.cc), **{k: np.asarray(v) for k, v in params.items()})
        # rank[row] = position of display row `row` in the Rastermap order
        rank = np.empty_like(isort)
        rank[isort] = np.arange(isort.size)
        json.dump(dict(description="Rastermap ordering of viewer display rows. isort[k] is the "
                       "display row shown at position k; rank[row] is the inverse. cluster_id[row] "
                       "is the ALF cluster index.", params=params, timing=rep,
                       isort=isort.tolist(), rank=rank.tolist(),
                       embedding=[round(float(e), 4) for e in emb], cluster_id=cl.tolist()),
                  open(HERE / f"rastermap_{a.tag}.json", "w"))


if __name__ == "__main__":
    main()
