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
  Within a tie the order is arbitrary: Rastermap sorts with numpy's default argsort, which
  does not keep input order for equal values. (An earlier version of this file said ties keep
  the input order; that was wrong.)
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

**Alternatives tried** (`alternatives.py`, `alternatives_hybrid.py`, results in
`alternatives.json`). All fitted on the first half and scored on the second half. S is the
neuron-level Rastermap-style score. "Within b" is the mean pairwise correlation inside blocks of
b adjacent neurons, which is what matters when a viewer cell averages b rows.

| Method | S | Within 3 | Within 10 | Within 30 | Within 100 | Seconds |
|---|---|---|---|---|---|---|
| Random | 0.0154 | 0.015 | 0.016 | 0.016 | 0.015 | 0 |
| Probe, shank, depth | 0.0189 | 0.069 | 0.063 | 0.051 | 0.033 | 0 |
| Rastermap, seeds 0 to 2 | 0.0259 to 0.0268 | 0.094 to 0.096 | 0.089 to 0.091 | 0.081 | 0.063 to 0.064 | 24 to 30 |
| Rastermap, 150 clusters | 0.0265 | 0.101 | 0.094 | 0.083 | 0.062 | 109 |
| Rastermap, n_splits = 1 | 0.0259 | 0.101 | 0.094 | 0.080 | 0.062 | 34 |
| Rastermap + tie-break | 0.0262 | 0.100 | 0.094 | 0.082 | 0.065 | +0.2 |
| Rastermap + refinement | 0.0262 | 0.140 | 0.112 | 0.089 | 0.065 | +2.4 |
| Sort by first PC | 0.0261 | 0.046 | 0.046 | 0.047 | 0.044 | 7 (PCs) |
| First PC + refinement | 0.0256 | 0.137 | 0.110 | 0.089 | 0.068 | 7 + 10 |
| Spectral ordering | 0.0260 | 0.043 | 0.043 | 0.042 | 0.041 | 7 + 7 |
| Average-linkage tree | 0.0243 | 0.149 | 0.113 | 0.085 | 0.063 | 7 + 1.3 |
| Tree + optimal leaf order | 0.0242 | 0.155 | 0.116 | 0.088 | 0.063 | 7 + 317 |
| Tree, flipped by Rastermap position | 0.0253 | 0.151 | 0.114 | 0.087 | 0.064 | Rastermap + 2 |

- "Refinement" is my own pass, not part of Rastermap: smooth the neurons' PC loadings along the
  current order, move every neuron to the position whose smoothed profile it matches best, and
  repeat with a shrinking smoothing width (25, 12, 6, 3 neurons).
- "Flipped by" keeps the tree and, at every branch, puts first the side with the lower mean
  Rastermap position.
- S barely separates methods: sorting by the first PC alone scores as high as Rastermap but has
  half its local correlation. S is dominated by the large-scale layout.
- Rastermap's weak spot is the local order (it only optimizes the order of 100 clusters of about
  43 neurons). Trees are best locally and somewhat worse in large-scale layout.
- Caveats: one session, one split, no error bars beyond the three Rastermap seeds. Within-3
  values could be raised by duplicate or split units, which trees would place side by side.
- Rastermap's time is mostly the cluster ordering: 15 of 24 s at 100 clusters, 1.2 s at 50,
  96 s at 150.

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

# Later on 2026-11-22: the lab's fork, new orderings, search, and slow drift

The lab now has a fork of Rastermap, `D:\Dropbox\code\rastermap`
(github.com/SteinmetzLab/rastermap, branch `fast-sort-and-local-order`; see `STEINMETZLAB.md`
there). It holds the refinement and tree orderings (`postsort.py`) and a faster cluster-order
search (`fastsort.py`), with tests. The scripts below use it (environment
`D:\temp\rastermap_env`, fork installed with `pip install -e`).

## Orderings shown in the viewer (`make_orderings.py` -> `orderings.json`)

Six whole-session orderings (100 ms bins, Rastermap defaults, seed 0):

| Key | What it is |
|---|---|
| `rm_rastermap` | Plain Rastermap. Identical to `rastermap_main.json` (asserted). |
| `rm_refine` | Rastermap, local order refined (`postsort.refine_isort`). |
| `rm_tree` | Average-linkage tree, branches flipped to follow Rastermap (`postsort.tree_isort`). |
| `hp_rastermap`, `hp_refine`, `hp_tree` | The same three, fitted after removing every cosine mode with period over 150 s from each neuron's binned counts (see "Slow drift"). |

`orderings_check.png` shows all six on the whole session (10 s bins, rate / session mean,
blocks of 5 neurons averaged). The viewer never filters the spikes; only the sorting input is
filtered.

## Can a better search find a much better global score? No. (`fast_search.py`, `fast_search.json`)

Fixed clusters (k-means seed 0, first half of the session), Rastermap's own objective F:

| Search | Wall time (16 threads) | F |
|---|---|---|
| Start (k-means order) | | 264.5 |
| Rastermap's greedy search | 14 to 15 s | 321.56 |
| Annealing, 1 run of 200,000 moves | 2.3 s | 325.1 |
| Annealing, 16 runs of 200,000 moves | 5 s | 322.3 to 325.7 (median 324.7) |
| Annealing, 16 runs + exhaustive polish | 29 s | 326.54 |
| 1,000 annealing runs of 200,000 moves | 490 s | 317.6 to 325.9 (median 324.6); top 5 polished: 326.26 to 326.54 |
| 160 annealing runs of 2,000,000 moves | 600 s | 325.2 to 326.3; top 5 polished: 326.53 to 326.55 |

- The greedy result is a true local optimum (an exhaustive sweep finds no improving move), so
  Rastermap does converge. 977 of 1,000 annealing runs beat it.
- The ceiling is about 326.55, 1.5% above the greedy result. The two best orders, from
  independent experiments, are the same layout (Spearman 0.999, 89 of 99 adjacent cluster pairs
  shared). The greedy order differs from them more (Spearman 0.86, 64 of 99 pairs).
- My earlier estimate that incremental evaluation would be about 100 times cheaper was wrong for
  full sweeps: a move that shifts m positions costs O(m n), and the average move in an exhaustive
  sweep shifts half the clusters. The saving is real only for the local moves annealing mostly
  makes. One exhaustive sweep costs about 1.3 s, so the polish is slow (about 1.5 s per move).
- These runs used a start temperature of 0.05; the fork's default is now 0.3, which was safer
  on unstructured test matrices and within 0.3% on this problem.

Across 100 k-means seeds (annealing, 16 runs each, no polish), scored at the neuron level on
the held-out half (S as defined above):

| | Held-out S |
|---|---|
| Annealing, 100 seeds | 0.0258 to 0.0271 (median 0.0265) |
| Greedy, 20 seeds | 0.0259 to 0.0268 (mean 0.0264) |
| Paired difference, 20 seeds | +0.00013 ± 0.00007 (s.e.m.) |

- No seed is much better: the best of 100 is 2% above the median. Seeds that score well on the
  fitted half tend to score well held out (r = 0.71), so the spread is partly real.
- Annealing does not improve the neuron-level order: within-3 correlation is the same and
  within-30 is 0.6% lower (0.0804 vs 0.0809, paired, 4 of 20 seeds better).
- Across seeds, greedy is usually close to annealing's cluster objective (ratio 1.004 on
  average, 0.996 to 1.036). Seed 0 is one of its worse cases.
- A better search does not make orders agree across seeds (`seed_consistency.py`): pairwise
  Spearman of neuron positions over 8 seeds is 0.75 (0.67 to 0.86) with annealing plus polish
  and 0.75 (0.67 to 0.83) with greedy. The differences between seeds come from the k-means
  clusters, not from the search. (Annealing plus polish was below greedy for 1 of the 8 seeds.)

## Slow drift (`drift.py`, `drift.json`)

Seen in the viewer in normalized mode: Rastermap groups neurons whose rate declines over the
session, rises over it, or peaks mid-session. Time courses were decomposed into cosine modes
(DCT-II; mode k has period 2T/k, and each mode is one number per neuron). Mode 1 is a monotone
drift, mode 2 a peak or dip in the middle.

Whole session (z-scored 100 ms counts, population mean projected out as Rastermap does). The
null keeps each neuron's own spectrum and flips the sign of each of its modes at random, i.e.
independent neurons with the same slow statistics:

| Periods | Modes | Share of variance in the top 10 PCs | SD of pairwise correlation | Same, null |
|---|---|---|---|---|
| over 2,400 s (modes 1, 2) | 2 | 26% | 0.769 | 0.769 |
| 1,214 to 1,619 s | 2 | 1.2% | 0.720 | 0.718 |
| 607 to 972 s | 4 | 0.4% | 0.517 | 0.504 |
| 304 to 540 s | 8 | 2.3% | 0.469 | 0.381 |
| 152 to 286 s | 16 | 2.7% | 0.384 | 0.260 |
| 76 to 147 s | 32 | 9.9% | 0.434 | 0.190 |
| 38 to 75 s | 64 | 5.8% | 0.295 | 0.132 |
| 19 to 38 s | 128 | 4.5% | 0.197 | 0.091 |
| under 1.2 s | 20,191 | 26% | 0.022 | 0.007 |

(Rows for 1.2 to 19 s are in `drift.json`.) Modes 1 and 2 hold 2.7% of a single neuron's
variance, 10% of the variance in the 200 PCs Rastermap sorts on, and 26% of the variance in its
top 10 PCs. Their pairwise correlations are spread exactly as for independent neurons.

Split halves: is a pair's within-band correlation in the first half reproduced in the second?
(Pearson r across the 9.3 M pairs)

| Periods | Modes per half | Reliability |
|---|---|---|
| over 1,214 s | 2 | 0.05 |
| 607 to 810 s | 2 | 0.05 |
| 304 to 486 s | 4 | 0.27 |
| 152 to 270 s | 8 | 0.37 |
| 76 to 143 s | 16 | 0.67 |
| 38 to 74 s | 32 | 0.68 |
| 0.6 to 38 s (six bands) | 64 to 2,048 | 0.49 to 0.62 |
| under 0.6 s | 8,047 | 0.69 |

So correlation structure is as reproducible as it gets once a band holds about 16 modes or more,
and not at all with 2. For the whole session, 32 modes in the lowest octave kept means a
cutoff period near 150 s, which is what the `hp_` orderings use.

Effect of the cutoff (fit on the first half after high-passing, anneal sorting; scored on the
second half):
- Held-out local scores barely depend on the cutoff. Within-3 correlation of Rastermap is 0.096
  to 0.099 for no filter and for cutoffs of 1,000, 300, 100, 30 and 10 s. Scored on fast
  activity only (held-out half high-passed at 30 s), it is 0.084 with no filter and 0.084 to
  0.089 with one (highest at the 10 s cutoff). So drift mostly changes the large-scale layout,
  not who sits next to whom.
- The refinement and tree gains hold with drift removed: at a 100 s cutoff, scored on the 30 s
  high-passed half, within-3 is 0.085 (Rastermap), 0.129 (refined), 0.140 (tree).
- Slow drift is not purely chance. Neighbors in an order fitted only on periods under 10 s in
  the first half still share slow components (periods over 300 s) in the second half: mean
  within-30 correlation of the slow parts 0.125, against 0.013 for a random order and 0.077 for
  depth order. An order fitted on unfiltered data gives 0.17.

## What `time_lag_window` does (from the source, `cluster.compute_cc_tdelay`)

It changes only the similarity between clusters used to order them. With `time_lag_window = w`
(in time bins), the similarity of cluster i to cluster j is the largest correlation between the
mean activity of i shifted by 0, 1, ..., w bins and that of j (i later than j), not the
correlation at zero lag. With `symmetric=False` (default) the matrix is asymmetric, and since
the objective uses only pairs (earlier position, later position), clusters whose activity
follows one another get placed next to each other in a consistent temporal direction. In the
code as written, the cluster that fires later gets the lower position. The k-means step and the
placement of single neurons still use zero-lag similarity. Not used in any ordering here.
