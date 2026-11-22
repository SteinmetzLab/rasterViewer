# Multiscale population-raster viewer for FD_016 2026-11-13/1: plan

Written 2026-11-21. Data: `D:\temp\fd016_share\alf` (4,325 QC units, 7 NP2.0 Quad probes,
60,997,315 spikes, 1.97 to 2430.81 s on the Timeline clock). Context on the data: the movie
README at `D:\Dropbox\presentations\2026-11-20_BreakingBarriers\rasterMovie\README.md`.

## The request

A viewer of the neurons x time matrix that opens on the whole session and lets you zoom in
"Google Maps style". Each screen pixel shows the mean firing rate over every neuron and every
instant it covers. At deep zoom this becomes a raster, because pixels hold 0 or 1 spikes.

## Verdict: yes, it makes sense and is very doable

The data are small. All 61 M spikes, stored exactly at 30 kHz resolution, compress to
about 115 MB. A complete precomputed count pyramid adds about 124 MB (both measured, see
"Data format"). The whole thing fits in a browser tab and in GPU memory, so neither version
needs tile streaming or a server. The hard parts are display choices, not engineering.
Those choices are below.

## What the example stills show (`figs/`)

All stills come from `render_view.py`. Each image pixel is one screen pixel, and its value is
the exact mean rate over the neurons and time it covers. No intermediate binning:
each spike goes to a pixel column by its own time. Neurons are weighted by how much of their
row falls inside the pixel. `test_render_view.py` checks this against a brute-force per-pixel
loop on 6 random views, and checks spike conservation on the full view (all 60,997,315
spikes accounted for). All pass. Rows are ordered by probe, shank, then depth (shallow on top),
the same order the movie used. Images are 1400 x 800 pixels.

| File | View | 1 pixel covers | Pixels with 0 spikes |
|---|---|---|---|
| `01_full_session_all_neurons.png` | Whole session, all 4,325 neurons | 5.4 neurons x 1.73 s | 0.7% |
| `02_..._normalized.png` | Same, each neuron divided by its session mean rate | same | same |
| `03_60s_all_neurons.png` | 60 s, all neurons | 5.4 neurons x 43 ms | 35% |
| `04_4s_all_neurons.png` | 4 s, all neurons | 5.4 neurons x 2.9 ms | 91% |
| `05_full_session_one_shank.png` (+ `06_..._normalized`) | Whole session, 93 neurons (Probe 3 shank 0) | 0.12 neurons x 1.73 s | 19% |
| `07_4s_one_shank.png` | 4 s, 93 neurons | 0.12 neurons x 2.9 ms | 98% |
| `08_500ms_60_neurons.png` | 500 ms, 60 neurons | 0.075 neurons x 0.36 ms | 99.8% |
| `09_50ms_20_neurons.png` | 50 ms, 20 neurons, strict pixel average | 0.025 neurons x 36 µs | 99.97% |
| `10_50ms_20_neurons_ticks.png` | Same, raster mode (spikes drawn 5 pixels wide) | same | |

The 4 s views use t = 1102 s, the movie's least-synchronized window. `07` shows a nice
traveling sweep across depth at around 1.6 to 1.9 s into that window.

## Five things the stills taught me, which shape the design

1. **Zoom has to be anisotropic.** Time spans a factor of about 2.4 million, from 2,429 s
   down to 1 ms. Neurons span only about 2 to 5 neurons per pixel at full view, even on a
   tall monitor. So time and neurons zoom independently, like a genome browser or an
   oscilloscope. The default would be wheel for time, shift+wheel for neurons, plus drag and
   box-zoom. Map-style uniform zoom would not work. Only the time axis needs a pyramid.
   Neuron aggregation is at most about 5 rows per pixel, cheap enough to do live.
2. **The raster regime starts much earlier than you might expect.** A pixel's expected spike
   count is (neurons per pixel) x (mean rate) x (pixel duration). Mean rate is 5.8 spikes/s.
   With all neurons on screen, that count reaches about 1 at 30 ms per pixel, a window of
   about 45 s. Every all-neuron view shorter than about a minute is therefore a quantized,
   salt-and-pepper raster, not a smooth rate map (`03`, `04`). This matches your definition
   exactly. It just means the "rate" reading only holds above about a minute. You may want an
   optional smoothing control for the in-between zooms. That would depart from the strict
   pixel-average rule, so it would be opt-in.
3. **Differences between neurons dominate the overview.** Session rates run from 0.007 to
   100 spikes/s, median 3. In `01`, the whole-session view is mostly horizontal stripes, each
   neuron's own mean rate. Dividing each neuron by its session mean (`02`) shows instead the
   population-wide state changes as vertical bands, which are the interesting part. I'd make
   this a toggle and default to normalized at coarse zoom. Normalization doesn't matter in
   raster mode.
4. **A single fixed log color scale works across all zooms.** I used 0.3 to 300 spikes/s.
   Coarse views use the middle of the range, and raster pixels simply saturate to the top
   color. No auto-contrast was needed, and fixed is better for comparing views.
5. **Below about 1 ms per pixel, "rate" stops meaning anything.** One spike in a 36 µs pixel
   is "28,000 spikes/s", and it is 1 pixel wide and nearly invisible (`09`). Below roughly
   1 to 2 ms per pixel, the viewer should switch to raster mode: each spike drawn as a tick of
   at least N pixels at its exact time (`10`). A raster mode is useful anyway, so you can see
   exact spike timing down to 33 µs, the 30 kHz sampling.

Smaller points:
- "Average over neurons" depends on row order. Within a shank, adjacent rows are neighbors in
  depth, so averaging is meaningful. At shank and probe boundaries a pixel can mix two shanks.
  An option is to snap pixel rows so they never straddle a shank boundary.
- Total spikes / (neurons x time) equals the mean of the per-neuron rates, so "average rate
  over the neurons" and "pooled spikes" are the same thing. In normalized mode it is the mean
  of normalized rates.

## Data format (decided)

Two tiers, written once by a Python preprocessing script:

- **Exact spikes**, sorted by time: `ticks` (uint32, 30 kHz samples on the Timeline clock) and
  `row` (uint16, display row). Because they are sorted by time, any time window is one
  contiguous slice, found by binary search. Measured with Blosc-zstd and bitshuffle:
  delta-encoded ticks 18 MB, rows 96 MB, about 115 MB total (raw 366 MB).
- **Count pyramid**: neurons x time, uint16 counts, base bin 16 ms, halving the resolution in
  time only at each level, 10 levels (16 ms to 8.2 s). Measured: 124 MB compressed total
  (36 MB for the finest level).
- Plus a small table of rows: cluster id, probe, shank, depth, session mean rate.

Container: plain **Zarr v3** (chunked, compressed arrays, readable from Python with
`zarr-python` and from JavaScript with `zarrita.js`). On the options you mentioned:

- **OME-Zarr** adds image metadata on top of Zarr. It would let napari or Neuroglancer open
  the pyramid as an image. But those viewers pick the nearest pyramid level and sample it.
  They do not average over each pixel's footprint, so their zoomed-out views would alias
  (skip bins) and break your rule. They also can't use the spike tier. I'd optionally add the
  OME-Zarr `multiscales` attributes, since it's free and gives a napari sanity check, but not
  build on it.
- **Data lakes and table formats** (Parquet, DuckDB, Iceberg) are for querying many sessions,
  not for drawing pixels. They're not needed at this scale. They become relevant if this turns
  into a viewer across many sessions.
- **No tiling/streaming**: at 240 MB total, load everything at startup. The format stays
  chunked, so streaming can be added later for bigger datasets.

How the tiers are used: if the visible window holds fewer than about 3 to 5 M spikes, draw
from the exact spikes. With all neurons, that is windows under 2 to 3 minutes; with fewer
neurons it covers longer windows. Otherwise use the pyramid. Pick the level whose bin is
about 1/16 of the pixel duration or finer, and snap pixel edges to that level's bin edges, so
pyramid views are exact too. Zoom steps are then quantized by at most about 6%, which is not
noticeable.

## Plan A: browser (WebGL2)

- **Rendering, spike tier.** All spikes go to the GPU once (366 MB raw, fine on a desktop
  GPU, tight on a laptop's integrated graphics). Each frame draws the visible contiguous
  range as points, with additive blending into a float texture: one texel per (pixel column,
  neuron row). Each spike is counted once, at its exact time. A second pass applies the
  vertical row-overlap weights, divides by neurons x duration (or applies per-neuron
  normalization), and maps through the color scale. This is the same math as
  `render_view.py`. Drawing about 5 M points per frame is comfortably interactive.
- **Rendering, pyramid tier.** Levels are stored as integer textures. The 151,926-bin finest
  level exceeds the 16,384-texel texture width limit, so it is wrapped or split. A fragment
  shader sums each pixel's footprint, about 16 bins x up to 5 rows.
- **Raster mode** below about 1 to 2 ms per pixel: spikes drawn as ticks at least N pixels
  wide (instanced quads).
- **User interface**: independent time and neuron zoom, a whole-session overview strip with
  the current window marked, probe and shank bands on the left, hover readout (cluster id,
  probe, shank, depth, time, rate, spike count), normalization toggle, color scale limits,
  jump to time, and URL-encoded view state (so a view can be bookmarked or sent).
- **Validation**: read the rendered float buffer back out of the page in headless Chrome
  (which I can drive here) and compare pixel-for-pixel against `render_view.py` over many
  random views, both tiers. A naive-then-critical review of the shader math by subagents.
- **Hosting**: a local static server (`python -m http.server`) is enough for you. Netlify
  would work technically, but the URL is public. This is unpublished data, so I won't deploy
  it publicly without your say-so.
- **Effort**: roughly a day of my work for a solid version, plus your feedback rounds.

## Plan B: not restricted to a browser

- **B1, Python desktop app** (Qt with vispy or pyqtgraph). Everything sits in RAM (68 GB on
  this machine). Each redraw slices the visible spikes and bins them with numpy/numba (about
  20 to 50 ms for 5 M spikes), or reads the pyramid. It reuses the tested `render_view.py`
  directly, and it's easiest to hook into analysis: click a neuron to get its waveform,
  autocorrelogram and peri-stimulus histogram; overlay behavior; export exact figures. The
  cons: about 10 to 30 redraws per second on big views, not smooth 60, and to share it others
  must install a Python environment.
- **B2, Datashader + HoloViews/Panel** (Python server, browser front end). The least code,
  maybe a couple of hours, and Datashader's per-pixel aggregation is close to your rule. (It
  assigns a spike to a pixel row by the row's center rather than by area; I'd fix the
  denominator.) But each zoom re-renders on the server in about 0.2 to 1 s, so it feels
  sluggish. It's good for trying display choices quickly, not as the final tool.
- **B3, native GPU app** (C++ or Rust). Fastest, but unjustified at this data size.

## Tradeoffs

| | Browser WebGL (A) | Python desktop (B1) | Datashader (B2) |
|---|---|---|---|
| Smoothness | 60 fps | 10 to 30 fps | 1 to 5 redraws/s |
| Share with others | A link or a folder | Needs a Python environment | Needs a running server |
| Exact per-pixel averaging | Yes (validated against reference) | Yes (is the reference) | Nearly |
| Link to analysis code | Harder (JavaScript) | Easy | Easy |
| Build effort | About a day | Half a day | Hours |
| Scales to many sessions | Yes, with streaming from Zarr | Limited by RAM | Yes, slower |

**My recommendation:** Plan A. The data are small enough that the browser has no real
disadvantage, it is the smoothest option, and it's the only one you can hand to someone
else. If you want something to play with tonight, B2 first, then A.

## Open questions for you

1. Who is it for? Just you exploring, the lab, or a public or talk audience? This decides
   hosting, and whether unpublished data may go on a public URL.
2. Default coloring: absolute rate (log scale) or normalized by each neuron's mean? Or some
   other normalization (z-score, rate relative to a baseline)?
3. Behavior events (stimulus onsets, choices, wheel movements) would be very valuable as tick
   marks along the time axis. They aren't in the share folder. Do they exist for FD_016 on
   the Timeline clock?
4. Row order: probe, shank, depth (as now), or an option for a data-driven order like
   Rastermap (Stringer, Pachitariu and colleagues; I did not check the exact citation)? That
   ordering would make the whole-session view far more informative than depth order.
5. This sits in a folder called `sonify`. Should the viewer eventually drive a sonification
   (for example, playing the spikes in the current view)? That affects the design of the
   playback and cursor.
