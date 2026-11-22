# Rastermap ordering of FD_016 (computed 2026-11-22)

A one-time Rastermap sort of all 4,325 units of FD_016 2026-11-13/1. The viewer's "Rastermap
order" mode uses `rastermap_main.json`, copied to `../data/rastermap.json` by
`../prep_rastermap.py`.

## The ordering to use

`rastermap_main.json` (also `rastermap_main.npz`, same content):
- `isort[k]`: the viewer display row shown at position k (rows are the viewer's rows: probe,
  shank, then y along the shank).
- `rank[row]`: the inverse, the position of display row `row`.
- `embedding[row]`: Rastermap's continuous 1D coordinate.
- `cluster_id[row]`: ALF cluster index of each display row.
- `params`, `timing`: settings and run stats.

Settings: Rastermap 1.0 (MouseLand, `pip install rastermap`) with its defaults: n_clusters 100,
n_PCs 200, locality 0 (global sort), time_lag_window 0, grid_upsample 10, per-neuron z-scoring
and the population-mean projection removed (both internal to Rastermap), k-means seed 0.
Input: spike counts in 100 ms bins over the whole session (t = 1.9685 to 2430.7685 s,
24,288 bins), all 4,325 units, all 61.0 M spikes. Built from the viewer's spike data
(`../data`; the run used the raw ticks.u32 and rows.u16, which the packed files reproduce
exactly), so it needs nothing from the Windows desktop.

## Cost (measured on the desktop)

The fit takes about 20 s regardless of length (18.5 to 34 s from 2,400 to 97,000 time bins),
and peak memory is 2.7 GB at 100 ms bins (6.9 GB at 25 ms, mostly my binning step). It would run
fine on the laptop too. Subset timing runs were done first (5% to 50% of the session).

## Checks (`validate_rastermap.py`, results in `validation.json`)

**It does put similar neurons next to each other, on data it did not see.** Fitted on the first
half of the session and scored on the second half, the mean correlation (100 ms bins) between a
neuron and its neighbor k positions away:

| Order | k = 1 | k = 5 | k = 10 | k = 50 |
|---|---|---|---|---|
| Rastermap, fit on first half | 0.097 | 0.088 | 0.080 | 0.049 |
| Probe, shank, depth (viewer) | 0.071 | 0.059 | 0.048 | 0.022 |
| Random | 0.016 | 0.016 | 0.016 | 0.015 |

Absolute correlations are small because single neurons at 100 ms are noisy. Depth order is
already good at k = 1 (neighbors on a shank share local activity, and possibly some split or
contaminated units); Rastermap's advantage grows at larger k, which is the regime that matters
when a pixel row aggregates many neurons.

**The fine-scale order is not unique.** Compared with the main ordering:

| Variant | Spearman of positions | Share of 20 nearest neighbors kept |
|---|---|---|
| k-means seed 1 | 0.60 | 16% |
| 250 ms bins | 0.56 | 14% |
| 25 ms bins | 0.41 | 11% |
| Fit on the first half only | 0.65 | 13% |

Chance for the neighbor share is 0.5%, so the variants agree far above chance, but the
large-scale arrangement of the 100 clusters and the order within them shift between runs.
Treat the main ordering as one good ordering, not the ordering. The variants' files
(`rastermap_seed1`, `_bin25ms`, `_bin250ms`, `_firsthalf`) are kept for comparison.

**Why orderings differ, and whether it matters** (`uniqueness.py`, results in
`uniqueness.json`). Rastermap does not optimize the neuron order directly. It runs three stages:
1. k-means of the neurons into 100 clusters, which is non-convex and starts from the seed;
2. a greedy search over the order of the 100 clusters, maximizing sum cc_ij x (-log distance);
3. placement of each neuron at the best-matching point on a 1,000-point grid between the
   ordered clusters.

Findings:
- Stage 2 stops at local optima. From 30 random starting orders on the same clusters, the
  objective ended between 316.5 and 326.5, and 20 of 30 beat the default run (321.6). The cluster
  orders found still agree with each other (Spearman 0.74 to 0.93).
- Stage 1 gives a different problem for each seed, so stage 2's objective can't be compared
  across seeds.
- Stage 3 ties: only about 740 distinct positions for 4,323 neurons (about 5.8 per position).
  Within a tie, the order is just the input order (probe, shank, y).
- Scores: I scored 10 seeds (fit on the first half) with a neuron-level version of the same
  objective: the mean over neuron pairs of correlation x the -log distance weight, normalized to
  mean 1.
  - On the held-out second half, seeds range from 0.0259 to 0.0268. Depth order scores 0.0189
    and random orders 0.0154.
  - So the spread between seeds is about 8% of Rastermap's gain over random, and about 12% of its
    gain over depth order.
  - Seed rankings on the two halves broadly agree (top four 8, 5, 4, 3 on the training half;
    1, 8, 4, 5 held out). The small differences are therefore partly real.
  - Orders fitted to the same data agree across seeds at Spearman 0.66 to 0.83.
- Conclusion: the score landscape is flat. Orderings that are quite different neuron by neuron
  score within a few percent of each other.

**Composition.** 48% of adjacent pairs come from different probes (0.1% in depth order), so a
Rastermap mode in the viewer would have to drop the rule that cells never aggregate across shank
or probe boundaries. One block of the order (positions about 2,600 to 3,000) is almost all probes
0 and 1. Position correlates with log session rate (Spearman -0.31): low-rate units, which are
noisiest after z-scoring, tend toward one end.

`rastermap_check.png`: the whole session at 1 s bins (z-scored per neuron), depth order above and
Rastermap order below. The colored strip at the left marks probe (tab10 colors, probe 0 blue,
1 orange, 2 green, 3 red, 4 purple, 5 brown, 6 pink).

## Files

| File | Role |
|---|---|
| `run_rastermap.py` | Bins spikes, fits Rastermap, writes `rastermap_<tag>.{json,npz}` |
| `validate_rastermap.py` | Stability, held-out neighbor correlation, composition, figure |

Environment: `D:\temp\rastermap_env` (uv venv, Python 3.12, rastermap 1.0, numba 0.68,
numpy 2.5, matplotlib). On the laptop, `pip install rastermap psutil matplotlib` is enough.
