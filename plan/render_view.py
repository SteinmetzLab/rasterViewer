"""Reference renderer for a multiscale population raster of FD_016 2026-11-13/1.

Given a view (time span x neuron-row span) and a pixel grid, each pixel's value is the mean
firing rate over every neuron and every instant it covers:

    rate(pixel) = sum_i w_i * count_i(pixel time span) / (sum_i w_i * dt_pixel)

where w_i is the fraction of neuron row i that falls inside the pixel's vertical extent
(rows are unit-height intervals [i, i+1) in display order). Time is exact: each spike
lands in exactly one pixel column by its own time, with no intermediate binning. This is
the oracle a viewer's output should match.

Display order: probe, shank, then depth descending (median Kilosort spike_positions y from
the movie export, D:/temp/fd016_raster/units.npz, whose unit index equals the ALF cluster
index; checked: per-unit counts in 1102-1182 s agree exactly).
"""
from pathlib import Path
import numpy as np

ALF = Path(r"D:/temp/fd016_share/alf")
MOVIE_NPZ = Path(r"D:/temp/fd016_raster/units.npz")
T_START, T_END = 1.9685, 2430.811   # synced period, seconds (README)


class Raster:
    def __init__(self):
        st = np.load(ALF / "spikes.times.npy")
        sc = np.load(ALF / "spikes.clusters.npy")
        meta = np.load(MOVIE_NPZ)["meta"]           # probe, shank, median y (um)
        self.n = meta.shape[0]
        self.order = np.lexsort((-meta[:, 2], meta[:, 1], meta[:, 0]))  # row -> cluster
        self.row_of = np.empty(self.n, np.int64)
        self.row_of[self.order] = np.arange(self.n)
        self.probe_of_row = meta[self.order, 0].astype(int)
        self.shank_of_row = meta[self.order, 1].astype(int)
        self.t = st                                   # sorted ascending
        self.r = self.row_of[sc].astype(np.int32)
        dur = T_END - T_START
        self.mean_rate_row = np.bincount(self.r, minlength=self.n) / dur

    def view(self, t0, t1, n0, n1, W, H, normalize=False):
        """Return (H, W) mean rate in spikes/s (or rate / each neuron's session mean)."""
        i0, i1 = np.searchsorted(self.t, [t0, t1])
        t = self.t[i0:i1]; r = self.r[i0:i1]
        r_lo, r_hi = int(np.floor(n0)), int(np.ceil(n1))
        keep = (r >= r_lo) & (r < r_hi)
        t = t[keep]; r = r[keep] - r_lo
        nr = r_hi - r_lo
        col = np.minimum(((t - t0) / (t1 - t0) * W).astype(np.int64), W - 1)
        if normalize:
            wts = 1.0 / np.maximum(self.mean_rate_row[r + r_lo], 1e-9)
        else:
            wts = None
        counts = np.bincount(r * W + col, weights=wts, minlength=nr * W).reshape(nr, W)
        # vertical area weights: overlap of row interval [i, i+1) with pixel [y_k, y_k+1)
        edges = np.linspace(n0, n1, H + 1)
        lo = np.maximum(edges[:-1, None], np.arange(r_lo, r_hi)[None, :])
        hi = np.minimum(edges[1:, None], np.arange(r_lo, r_hi)[None, :] + 1)
        Wv = np.clip(hi - lo, 0, None)                # (H, nr)
        dt = (t1 - t0) / W
        return (Wv @ counts) / (Wv.sum(1, keepdims=True) * dt)
