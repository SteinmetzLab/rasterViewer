"""Reference (CPU, float64) implementation of the viewer's pixel values.

The viewer's GPU output is tested against this. Semantics, for a grid of W x H cells
(a cell is block x block screen pixels):

Time. The view starts at T0 (in 30 kHz ticks, float) and each cell column spans dT ticks.
Column j's first integer tick is B_j = ceil(T0 + j*dT) (float64). A spike at integer tick s
is in column j iff B_j <= s < B_(j+1), i.e. column floor((s - T0) / dT) up to float64
rounding exactly at an edge. The viewer uses the same integer boundaries on the GPU.

Vertical axis (the "layout"). Every neuron occupies an interval [lo, hi) of a vertical
coordinate Y, and the Y axis is divided into segments. Three layouts:
- "rows" (sorted by depth): the neuron at display row i occupies [i, i+1). Rows are ordered by
  probe, shank, then y along the shank (shallow at top). Each (probe, shank) is a segment.
- "depth" (plotted by depth): Y is in micrometers, to scale. Each (probe, shank) is a segment,
  and every segment spans the same y range: from ymin - h/2 to ymax + h/2, where [ymin, ymax]
  is the recorded site range DEPTH_SITES_UM (0 to 2865 um, NP2.0 bank 0), widened if any unit
  lies outside it, and h = 15 um (the NP2.0 site pitch). Segments are stacked in probe, shank
  order with no gaps. A unit at height y above the shank's deepest site occupies a band of
  height h centered on y, shallow at top: within its segment, [top - y - h/2, top - y + h/2)
  with top = ymax + h/2. Units at similar depths overlap.
- "rastermap": the neuron shown at Rastermap position k (row isort[k]) occupies [k, k+1). The
  whole population is one segment, so cells aggregate across probes and shanks.
Cell row k (k = 0 at top) covers [n0 + k*dn, n0 + (k+1)*dn), dn = (n1 - n0) / H. The cell
belongs to the segment that contains its center, and its span is cut to that segment (and to
n1), so a cell straddling a segment boundary ignores the other segment. Each neuron is weighted
by the length of its interval's overlap with the cut span, divided by the interval's height
(1 for rows and rastermap, h for depth), i.e. by the fraction of the neuron inside the cell.

Value. Absolute mode: sum_i ov_i * c_i / (sum_i ov_i * dt_j), with dt_j = (B_(j+1) - B_j) / fs
the column's duration in whole samples (not dT / fs, which would bias columns that hold one
sample more or fewer than average), the mean rate in spikes/s over
the cell's neurons and time. Normalized mode: each neuron's count is first multiplied by
1 / (r0_i + eps), where r0_i is that neuron's session mean rate and eps a softening rate in
spikes/s, giving the mean over neurons of rate / (session mean + eps), dimensionless.
Special values: -1 if the column holds no 30 kHz samples; -2 if no neuron overlaps the cell
(only possible in the depth layout, between units).
"""
from pathlib import Path
import json
import numpy as np
from pack_data import load_spikes

DATA = Path(__file__).resolve().parent / "data"
DEPTH_H = 15.0      # um, height of a unit's band in the depth layout (NP2.0 site pitch)
DEPTH_SITES_UM = (0.0, 2865.0)   # deepest and shallowest site centers (NP2.0 bank 0, 192 x 15 um)


class Data:
    def __init__(self, data_dir=DATA):
        self.ticks, self.rows = load_spikes(data_dir)     # raw files, or the packed copies
        self.meta = json.load(open(data_dir / "meta.json"))
        r = self.meta["rows"]
        self.n = self.meta["n_rows"]
        self.fs = self.meta["fs"]
        self.r0 = np.array(r["r0"]); self.seg_start = np.array(r["seg_start"])
        self.seg_end = np.array(r["seg_end"])
        rm = data_dir / "rastermap.json"
        self.isort = np.array(json.load(open(rm))["isort"]) if rm.exists() else None
        self._layouts = {}

    def layout(self, kind, h=DEPTH_H):
        """(lo, hi, unit, order, seg_lo, seg_hi) for a layout. Position p shows row order[p]
        and occupies [lo[p], hi[p]); unit = hi - lo; segments are [seg_lo[s], seg_hi[s])."""
        key = (kind, h)
        if key in self._layouts:
            return self._layouts[key]
        n, segs = self.n, self.meta["segments"]
        if kind == "rows":
            lo = np.arange(n, dtype=np.float64); order = np.arange(n)
            seg_lo = np.array([s["start"] for s in segs], float)
            seg_hi = np.array([s["end"] for s in segs], float)
            unit = 1.0
        elif kind == "rastermap":
            lo = np.arange(n, dtype=np.float64); order = self.isort
            seg_lo, seg_hi = np.array([0.0]), np.array([float(n)])
            unit = 1.0
        elif kind == "depth":
            # integer units of 0.1 um, exactly as the viewer computes them (y is stored to 0.1 um)
            y10 = np.round(np.array(self.meta["rows"]["y_um"]) * 10).astype(np.int64)
            h10 = int(round(h * 10)); assert h10 % 2 == 0
            lo10 = np.zeros(n, np.int64); seg_lo10 = []; seg_hi10 = []; base = 0
            mn = min(int(round(DEPTH_SITES_UM[0] * 10)), int(y10.min()))
            mx = max(int(round(DEPTH_SITES_UM[1] * 10)), int(y10.max()))
            top10 = mx + h10 // 2; span10 = mx - mn + h10          # the same for every shank
            for s in segs:
                a, b = s["start"], s["end"]
                lo10[a:b] = base + top10 - y10[a:b] - h10 // 2
                seg_lo10.append(base); seg_hi10.append(base + span10); base += span10
            lo = lo10 / 10.0; order = np.arange(n); unit = h
            seg_lo, seg_hi = np.array(seg_lo10) / 10.0, np.array(seg_hi10) / 10.0
        else:
            raise ValueError(kind)
        out = (lo, lo + unit, unit, np.asarray(order), seg_lo, seg_hi)
        self._layouts[key] = out
        return out

    def extent(self, kind, h=DEPTH_H):
        return float(self.layout(kind, h)[5][-1])

    def cell_values(self, T0, dT, n0, n1, W, H, normalized=False, eps=0.1, layout="rows",
                    h=DEPTH_H):
        lo, hi, unit, order, seg_lo, seg_hi = self.layout(layout, h)
        Tend = T0 + W * dT
        i0 = np.searchsorted(self.ticks, np.ceil(T0), "left")
        i1 = np.searchsorted(self.ticks, np.ceil(Tend), "left")
        s = self.ticks[i0:i1].astype(np.float64); r = self.rows[i0:i1].astype(np.int64)
        bounds = np.ceil(T0 + np.arange(W + 1) * dT)     # first integer tick of each column
        col = np.clip(np.searchsorted(bounds, s, side="right") - 1, 0, W - 1)
        counts = np.bincount(r * W + col, minlength=self.n * W).reshape(self.n, W).astype(np.float64)
        if normalized:
            counts *= (1.0 / (self.r0 + eps))[:, None]
        counts_pos = counts[order]                         # rows in display (position) order
        dn = (n1 - n0) / H
        y0 = n0 + np.arange(H) * dn; y1 = y0 + dn
        yc = 0.5 * (y0 + y1)
        sg = np.clip(np.searchsorted(seg_lo, yc, side="right") - 1, 0, len(seg_lo) - 1)
        a = np.maximum(y0, seg_lo[sg]); b = np.minimum(np.minimum(y1, seg_hi[sg]), n1)
        ov = np.clip(np.minimum(b[:, None], hi[None, :]) - np.maximum(a[:, None], lo[None, :]),
                     0, None) / unit
        # each column holds (bounds[j+1] - bounds[j]) whole 30 kHz samples; a spike at tick s
        # stands for true times [s - 0.5, s + 0.5) / fs, so that is the column's duration
        dt_s = np.diff(bounds) / self.fs
        den = ov.sum(1, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = (ov @ counts_pos) / (den * dt_s[None, :])
        v = np.where(den > 0, v, -2.0)
        return np.where(dt_s[None, :] > 0, v, -1.0)     # -1 = column holds no samples


if __name__ == "__main__":
    # self-check against a slow per-cell loop (rows layout, the original definition)
    D = Data()
    rng = np.random.default_rng(1)
    for trial in range(5):
        W, H = int(rng.integers(10, 40)), int(rng.integers(5, 30))
        T0 = rng.uniform(1e5, 7e7); dT = rng.choice([1.5, 30.0, 3333.3, 1e5])
        n0 = rng.uniform(0, 4000); n1 = min(D.n, n0 + rng.choice([6.5, 60.2, 900.0]))
        norm = bool(trial % 2)
        v = D.cell_values(T0, dT, n0, n1, W, H, norm, 0.1)
        ref = np.zeros((H, W)); dn = (n1 - n0) / H
        m = (D.ticks >= T0) & (D.ticks < T0 + W * dT)
        s = D.ticks[m].astype(float); r = D.rows[m].astype(int)
        for k in range(H):
            y0, y1 = n0 + k * dn, n0 + (k + 1) * dn
            c = int(np.floor((y0 + y1) / 2)); a, b = max(y0, D.seg_start[c]), min(y1, D.seg_end[c])
            for j in range(W):
                mm = (s >= T0 + j * dT) & (s < T0 + (j + 1) * dT)
                num = den = 0.0
                for i in range(int(np.floor(a)), int(np.ceil(b))):
                    ov = min(b, i + 1) - max(a, i)
                    if ov <= 0: continue
                    w = 1 / (D.r0[i] + 0.1) if norm else 1.0
                    num += ov * w * np.sum(r[mm] == i); den += ov
                nt = np.ceil(T0 + (j + 1) * dT) - np.ceil(T0 + j * dT)
                ref[k, j] = num / (den * nt / D.fs)
        print("rows", trial, "max abs diff", np.abs(v - ref).max(), "max", ref.max())

    # depth and rastermap layouts against a slow loop written from the docstring's definition
    for kind in ["depth", "rastermap"]:
        lo, hi, unit, order, seg_lo, seg_hi = D.layout(kind)
        Y = seg_hi[-1]
        for trial in range(4):
            W, H = int(rng.integers(8, 25)), int(rng.integers(5, 25))
            T0 = rng.uniform(1e5, 7e7); dT = rng.choice([30.0, 3333.3, 1e5])
            span = Y * rng.choice([0.002, 0.03, 0.3, 1.0]); n0 = rng.uniform(0, Y - span)
            n1 = n0 + span; norm = bool(trial % 2)
            v = D.cell_values(T0, dT, n0, n1, W, H, norm, 0.1, layout=kind)
            ref = np.zeros((H, W)); dn = (n1 - n0) / H
            m = (D.ticks >= np.ceil(T0)) & (D.ticks < np.ceil(T0 + W * dT))
            s = D.ticks[m].astype(float); r = D.rows[m].astype(int)
            for k in range(H):
                y0, y1 = n0 + k * dn, n0 + (k + 1) * dn
                c = 0.5 * (y0 + y1)
                si = max(i for i in range(len(seg_lo)) if seg_lo[i] <= c or i == 0)
                a, b = max(y0, seg_lo[si]), min(y1, seg_hi[si], n1)
                for j in range(W):
                    b0, b1 = np.ceil(T0 + j * dT), np.ceil(T0 + (j + 1) * dT)
                    mm = (s >= b0) & (s < b1)
                    num = den = 0.0
                    for p in range(D.n):
                        ov = (min(b, hi[p]) - max(a, lo[p])) / unit
                        if ov <= 0: continue
                        row = order[p]
                        w = 1 / (D.r0[row] + 0.1) if norm else 1.0
                        num += ov * w * np.sum(r[mm] == row); den += ov
                    ref[k, j] = -2.0 if den == 0 else num / (den * (b1 - b0) / D.fs)
            print(kind, trial, "max abs diff", np.abs(v - ref).max(), "max", ref.max(),
                  "cells with no neurons", int((ref == -2).sum()))
