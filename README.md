# Population raster viewer

A zoomable view of all 4,325 neurons x 2,429 s of FD_016 2026-11-13/1 (7 NP2.0 Quad-Base
probes, Steinmetz Lab). It opens on the whole session. Every displayed cell is the exact mean
firing rate over every neuron and every instant it covers, so it turns into a raster by itself
as you zoom in. Online at https://www.steinmetzlab.net/rasterViewer/ (GitHub Pages, from this
repository's `main` branch). The design discussion and example stills are in
`plan/PLAN.md` and `plan/figs/`.

## Running it

Online: open the address above. Locally:

```
python serve.py       # then open http://localhost:8765/ (Windows: Chrome or Edge)
python3 serve.py      # the same on macOS (Chrome or Safari)
```

`serve.py` needs only the Python standard library (3.7 or later). Opening `index.html` as a
file does not work, because browsers block loading the data from `file://` pages.

Startup downloads 115 MB of compressed spikes (366 MB decoded; see "Data files"). Locally
that takes 1 to 2 s; online it depends on the connection (about 10 s at 100 Mbit/s).
It needs WebGL2 with `EXT_color_buffer_float` and `EXT_float_blend`, which desktop Chrome, Edge,
Firefox and Safari 17 have. If they are missing, the viewer refuses to start. A 16-bit fallback
would silently cap counts at 2,048. Tested in Chrome on the TITAN V (Windows), and in Chrome 158
and Safari 17.6 on an Apple M3 Max MacBook Pro (macOS 14).

**Mac laptop on battery:** Chrome's Energy Saver (Settings > Performance) caps every page at
30 frames/s when it is on. On the M3 Max MacBook it is set to turn on whenever the laptop is
unplugged. Plug in, or turn it off, for 60 to 120 frames/s. Safari ran at 60 frames/s on
battery.

## Using it

| Action | How |
|---|---|
| Zoom time | Mouse wheel over the raster or the time axis; on a trackpad, two-finger vertical scroll or pinch |
| Zoom neurons | Shift + wheel, or wheel over the neuron axis on the left |
| Pan | Drag, the arrow keys, or horizontal scrolling (sideways two-finger swipe, tilt wheel) |
| Box zoom | Ctrl + drag, Cmd + drag on a Mac, or right-drag |
| Jump in time | Click or drag in the population-rate strip at the top |
| Zoom to a probe or shank | Click its label (or its color band) on the neuron axis |
| Reset / all neurons / whole session | Buttons, or keys R / A / W; + and − zoom time |
| Sonify (play the view as sound) | The Sonify button, or S or the space bar; see "Sonify" below |
| Readout | Off by default. Click a cell to show it (it stays until the view changes; Esc clears it), or hold Alt (Option on a Mac) to have it follow the pointer. It gives the time, the neuron under the pointer (cluster id, probe, shank, position along the shank, session rate), and the cell's spike count, neuron count, duration and value |

Controls in the top bar:
- **Color by**: firing rate in spikes/s (default), or rate / (session mean + ε), where ε is a
  softening rate in spikes/s (default 0.1) set in the box beside it. Each mode keeps its own
  color range and map. The defaults are 0.3 to 50 spikes/s with magma, and 0.1 to 10 with
  viridis (a one-directional map, so 1, the neuron's own mean, is not visually special).
- **Neurons**: how neurons are placed on the vertical axis.
  - *Sorted by depth* (default): one row per neuron, ordered by probe, shank, then position
    along the shank, shallow at top.
  - *Plotted by depth (to scale)*: each shank is drawn to scale in micrometers, and each neuron
    is a 15 µm band (the NP2.0 site pitch) centered on its position, shallow at top. Neurons at
    similar depths overlap and are averaged; stretches with no neuron are drawn as a rate of zero. Every
    shank spans the same full recorded range (sites at 0 to 2865 µm above the deepest site,
    NP2.0 bank 0), so all shanks are drawn identically; depth below the brain surface still
    needs histology.
  - *Rastermap order*: the order in `../rastermap/rastermap_main.json` (see the README there),
    so neurons with similar activity are adjacent. Cells average across probes and shanks
    here. That ordering is one good ordering, not a unique one: reruns with another seed or
    bin size keep only about 15% of each neuron's 20 nearest neighbors. The colored strip on
    the left shows each neuron's probe.
  - Five further orderings (added 2026-11-22; `data/orderings.json`, built by
    `prep_orderings.py` from `rastermap/orderings.json`; details in `rastermap/README.md`).
    They are drawn exactly like Rastermap order, with another order:
    - *Rastermap, refined*: the Rastermap layout with the order of nearby neurons improved.
    - *Tree (Rastermap layout)*: a hierarchical-clustering tree of the neurons, its branches
      flipped to follow the Rastermap layout. Best local order of the three.
    - *Rastermap, no slow drift*, *Rastermap, refined, no slow drift*, *Tree, no slow drift*:
      the same three, computed after removing every component slower than a 150 s period from
      the sorting input. Slow drifts over the session then no longer drive the order. The
      spikes shown are never filtered.
- **Spikes**: which spikes count. *All* (default); *Tonic only*; *Burst spikes*; *Burst
  onsets* (only each burst's first spike, so the value is bursts/s); or *Burst + tonic*,
  which draws burst spikes black to magenta and tonic spikes black to green and adds the two,
  so a cell with plenty of both looks white. Both use the same color range; burst-spike rates
  are much lower than tonic rates, so bursts look dim when zoomed far out.
  A burst is a run of at least K spikes of one unit with every ISI at most M ms, whose first
  spike follows at least S ms without a spike from that unit. S, M and K are set in the boxes
  that appear beside the menu; the default (100 ms, 4 ms, 2) is the usual criterion for
  thalamic low-threshold bursts (Lu, Guido & Sherman 1992). Changing them re-runs detection
  in the browser (about 1 s). How this definition was chosen, with diagnostic figures, is in
  `bursts/README.md`.
- **Cell**: the aggregation unit, 1, 2 (default), 3, 4, 6 or 8 screen pixels square. Larger
  cells average over more neurons and time, and individual spikes become bigger blocks.
- **Color range** (log scale) and **Map**.
- The view is kept in the URL (`#t0=...&t1=...`), so a view can be bookmarked.

## Sonify

**Sonify** (the button, S, or the space bar) plays the view as sound. A line is drawn at the
middle time of the view, the raster scrolls left under it, and what crosses the line is played.
Press again to stop. Everything can be changed while it plays: pan, zoom, the neuron ordering,
the spike selection, the color range. Panning or zooming in time continues from the middle of
the new view.

- **Speed.** The raster always moves at the same speed on the screen: one plot width in 30 s
  (the *Crossing* setting: 10, 20, 30 or 60 s). The playback rate therefore follows the zoom:
  the whole 40 min session in view plays at 81 times real time, a 3 s window at 1/10 of real
  time. The status bar shows the rate.
- **Pitch.** The vertical axis of the current view maps to pitch, high at the top, over 4
  octaves centered on A4 (*Range*: 2 to 5). This holds for every neuron ordering and any
  vertical zoom. Each cell row sounds the nearest note of the chosen *Key* and *Scale* (major
  or minor pentatonic, major, natural minor, whole tone, 12-tone), or, with *No rounding*, one
  of up to 64 equal pitch bands. When neurons are taller than a cell, all of a neuron's cells
  sound the single note of the neuron's middle, so a neuron is one note.
- **Loudness follows color.** With t a cell's position on the color scale (0 at the bottom of
  the color range, 1 at the top), its level is 36 dB x (1 - t) below full, and silent at t = 0
  (no spikes, or a rate at or below the bottom of the color range). A note's level is the RMS
  of the levels of its cells. So the color range sets what is audible, and bright is loud.
- **Single spikes and firing rates are one mechanism.** Each note is an oscillator whose
  loudness follows the note's level with a fast attack and a slower release (the *Sound*
  setting: marimba, harp, organ, pure tone). Zoomed in, a spike is one bright cell: the level
  jumps for one column and falls back, which is a struck note that rings out. Zoomed out, the
  level changes a little from column to column: a sustained tone whose loudness follows the
  firing rate. In between, spikes that arrive faster than the release fuse into a tremolo.
- **What you hear is what you see.** The sound is computed from the cell values as drawn, one
  column at a time; the view advances by whole columns so a column holds the same spikes on
  every frame. Timing is therefore quantized to one column: 1/W of the crossing time, about 40
  ms of playing time for a 740-column plot at 30 s (finer with a wider plot or Cell = 1).
- Where the view cannot scroll (the window touches the start or the end of the session) the
  view stays and the line moves instead. Started with the window at the session start, the
  line starts at the left edge, so the beginning is not skipped.
- The settings are kept in the URL. Playing stops when the tab is hidden.

## What each cell shows (exact definition)

- **Time**: column j covers the integer 30 kHz ticks from ceil(T0 + j dT) up to
  ceil(T0 + (j+1) dT). Spike times were rounded to the nearest 30 kHz tick (±17 µs). They
  were already interpolated onto the Timeline clock, so sub-tick precision did not mean much.
- **Neurons**: each neuron occupies an interval of the vertical axis: [i, i+1) for row i when
  sorted by depth; a 15 µm band centered on its position when plotted by depth; [k, k+1) at
  Rastermap position k. A cell uses the neurons it overlaps, each weighted by the fraction of
  the neuron inside the cell, **but only from the shank under the cell's center** (in
  Rastermap order the whole population counts as one shank). A cell straddling a shank (or
  probe) boundary ignores the other shank's neurons. Within a shank, neurons are ordered by
  position along the shank (median Kilosort spike y, shallow at top), as in the Breaking
  Barriers movie. The readout shows that y as "µm above deepest site"; depth below the brain
  surface needs histology. The exact definition is in `reference.py`.
- **Value**: sum of overlap x count / (sum of overlap x column duration), i.e. the mean over the
  cell's neurons of each neuron's rate in that time. A column's duration is its whole number of
  30 kHz samples / 30,000, not the average window / W. When a column holds a fractional number
  of samples on average, columns alternate between n and n+1 samples, and dividing by the
  average would stripe the image. In normalized mode, each neuron's count is first divided by
  (its session mean rate + ε).
- **Shortest window**: W samples, so every column holds at least one sample. That is about
  47 ms at 1400 cells across; a larger Cell size allows shorter windows.
- A cell with no spikes gets the bottom color of the map. So does a cell with no neuron (only
  when plotted by depth); the readout under the cursor says which it is. A column with no
  30 kHz sample is mid gray.

## How it works

All spikes go to the GPU once, in two copies: sorted by time (so any window is one contiguous
range), and regrouped by neuron row, each row in time order. Each redraw has three passes:
1. Counts spikes into a texture with one texel per (cell column, neuron row), by one of two
   methods that give identical, exact counts. Both use integer column boundaries (ticks)
   computed in float64.
   - **points**: draws the window's spikes as 1-pixel points, additively blended. Cost grows
     with the number of spikes in the window.
   - **search**: each texel binary-searches its row's sorted spike ticks for the column's two
     boundaries; the count is the difference of the two positions. Cost grows with the grid
     (columns x rows), not with the spikes.

   Search was added for the Mac: Apple's tile-based GPUs are slow at rasterizing tens of
   millions of points. It is an optimization only, so the viewer protects against it going
   wrong:
   - **Startup check.** At every launch, both methods count the session's last 300 s for all
     neurons on a 512-column grid (about 0.2 s). The counts must be bitwise identical and add
     up to the spikes in the window. If they don't, or the shader fails to compile, or the GPU
     can't allocate the 244 MB texture, search is switched off and points is used for
     everything. The console says which.
   - **Automatic choice.** The same check times both methods. From the timings it sets the
     window density (spikes per texel) above which search is used. The M3 Max gave 0.13 to
     0.40 spikes per texel; the timings are noisy at this scale. Near the threshold both
     methods take a few ms, so the noise doesn't matter.
   - **Forcing a method.** Add `?count=points` or `?count=search` to the URL (before the `#`).
     `?count=points` also skips building the search data, which is the escape hatch if search
     ever misbehaves. The status bar shows which method drew the current view.
2. Averages each cell's rows within its shank, using row spans and fractional overlaps
   computed in float64 on the CPU.
3. Applies the log color map.

Burst modes: the browser classifies every spike once (one pass in time order), then keeps the
burst spikes and the burst onsets as two extra spike sets on the GPU (about 30 MB), each
countable by either method. Pass 1 then counts two sets (for example all spikes and burst
spikes), and pass 2 combines the two counts per texel before averaging, so tonic = all minus
burst is formed from exact integer counts. In the overlay, pass 2 outputs burst and tonic
values side by side and pass 3 adds their two color ramps.

The hover readout recomputes its cell on the CPU from the raw spikes.

Speed, whole session with all 61 M spikes and all neurons, GPU time per redraw:

| Machine | Grid | points | search |
|---|---|---|---|
| TITAN V, Chrome, Windows | 2488 x 1396 | about 42 ms | not yet measured |
| M3 Max, Chrome 158, macOS | 2488 x 1192 | 132 ms | 5 ms |
| M3 Max, Safari 17.6, macOS | 1992 x 1384 | 154 ms | 4 ms |

With automatic selection, the views timed on the M3 Max in Chrome (whole session down to
0.1 s, 20 to 4,325 neurons) took 4 ms or less of GPU time each. The whole frame, with axes and hover, is
about 10 ms. Before search, panning the whole session ran at about 6 frames/s on the Mac.

## Validation

- `reference.py`: float64 CPU implementation of the definition above. Its `__main__` checks it
  against a slow per-cell loop.
- `make_refs.py`: writes reference values for 45 views to `refs/` in this folder (about 185 MB,
  not in git: run it once before using `?test=1`; it reads the packed spikes). The first
  22 (all sorted by depth) cover:
  - whole-session views in both modes, a 3000 x 1500 grid, and a grid with Cell = 3;
  - cells straddling a shank boundary;
  - 50 ms and 30 ms deep zooms, including columns with no samples;
  - fractional column widths of 1.5 and 21 samples;
  - the session's last 2 s (largest tick values);
  - 10 random views.

  The other 17 cover plotted-by-depth views (whole session in both modes, one shank, a whole
  probe, a shank boundary, a 30 µm and a 2 µm deep zoom where cells are much thinner than a
  neuron's band) and Rastermap views (whole session in both modes, 200 neurons, 2 neurons at
  50 ms), plus 3 random views of each. The last 6 cover the burst modes with the default
  burst definition: burst spikes, tonic spikes (normalized), burst onsets on Probe 6 shank 0,
  the overlay in the window of the 1110.6 s burst episode (both of its outputs), burst spikes
  plotted by depth, and tonic spikes in Rastermap order. It also writes `refs/bursts.json`:
  the number of burst spikes and onsets, and the sums of their spike indices, from
  `bursts/detect_bursts.py`, which the self-test compares with the browser's own detection.
- `http://localhost:8765/?test=1` renders each view on the GPU, once with each counting
  method, and reads the values back. It then:
  - compares them with the references;
  - requires the two methods' raw counts to be bitwise identical, which is stricter than the
    1e-3 tolerance (an off-by-one count in a cell of 2,000 spikes is only 5e-4);
  - checks the CPU hover path at 40 random cells per view;
  - times the whole-session view with each method. The results are in the status bar, the console, and `window.__testResults`.
- Result on 2026-11-22 (Windows, TITAN V, points method only, which was all that existed):
  all 22 views have 0 cells differing by more than 1e-3 (relative). The largest relative error
  is 1.8e-5 (float32 accumulation). Total spikes match. The CPU hover path matches to 6e-8.
- Result on 2026-11-22 (M3 Max, Chrome 158 and Safari 17.6): the same for both methods: 0
  mismatched cells in all 22 views, largest relative error 1.8e-5, CPU hover path to 6e-8.
  The raw counts of the two methods are bitwise identical in all 22 views.
  **The search method has not yet run on the Windows machine.** Run `?test=1` there once. The
  startup check will switch search off there by itself if it misbehaves.
- Result on 2026-11-22 for the three vertical layouts (M3 Max, built-in Chromium): 0
  mismatched cells in all 39 views with both counting methods, raw counts identical between
  methods. Largest relative error 1.8e-5 sorted by depth, 7.5e-5 in Rastermap order. Plotted
  by depth, a few cells at exact band edges differ by about 3e-8 spikes/s: float64 rounds the
  edge in the reference and counts a 2e-12 µm sliver of the next neuron, while the viewer
  computes band positions in exact integer steps of 0.1 µm. Rewriting `reference.py` for the
  three layouts left the 22 original reference files unchanged.
- Result on 2026-11-22 with the burst modes (M3 Max, built-in Chromium, from the packed data):
  0 mismatched cells in all 45 views with both counting methods, raw counts identical between
  methods. The browser's burst detection matches the Python detector exactly (2,193,316 burst
  spikes in 787,718 bursts, and identical index sums). The self-test checks only the default
  definition; a one-off check of 50 ms, 6 ms, 3 spikes also matched exactly (2,297,895 burst
  spikes in 572,897 bursts), and an independent review compared the two detectors on 900
  random synthetic cases (K from 2 to 6, fractional thresholds, zero ISIs) with no mismatch.
- A second independent review (of the layouts and readout changes) found no value bug in
  normal cases. It found a rare float32 edge in the depth layout (a cell covered only by
  slivers of bands could show as empty; such cells now weight their slivers equally), and UI
  bugs, now fixed: about 5% of clicks failed to pin the readout because view clamping was
  not stable under repetition, a focused dropdown swallowed the keyboard shortcuts, and pin
  vs. Alt priority was inconsistent.
- An independent review of the Mac changes (a separate agent) found no value bug. It raised
  four risks, now handled:
  - a search failure on Windows could block startup;
  - a failed texture allocation could give silently wrong counts;
  - the switch threshold had been tuned on the Mac only;
  - the test tolerance could hide off-by-one counts.
- An independent critical review (a separate agent) found two value bugs, now fixed:
  - Dividing by the average column duration striped the image at fractional sample counts
    (2x at 1.5 samples per column). After the fix, at 1.5 samples per column, columns with 1
    and 2 samples average 6.55 vs 6.69 spikes/s (paired difference not significant, t = 1.3);
    before, 4.30 vs 8.74.
  - The R16F fallback capped counts at 2,048; it was removed.

  It also found smaller issues, now fixed: URL parameters were not validated, the last 2 spikes
  could not be shown, the colorbar was off by one color-table step, the raster could shift by a
  pixel at fractional display scaling, and "depth" was a misleading label.

- **Sonify**: `?sonifytest=1` runs 36 checks in the browser:
  - the notes of four scales against hand-written lists;
  - rows to notes for every ordering (highest note at the top, never rising downward, one note
    per neuron when zoomed in, and for "plotted by depth" against a brute-force search);
  - the transport (which column sounds when, and where the view and the line are) run on a
    simulated clock for six cases, including the session start, the session end, the whole
    session in view and stalled frames: every column is scheduled exactly once, before it
    sounds, from inside the view, and the line is over the column being heard;
  - the synthesis rendered offline with known levels (amplitude within 2% of level, attack and
    release times, silence before a note's column, no leakage between notes, frequencies);
  - for six real views, the note levels the player reads back from the picture against levels
    computed from the float64 readout (`cellCPU`) with separate code (largest difference 6e-8).

  Live runs in Chrome, recording the output with an analyser: the playback rate equals the
  window divided by the crossing time; the sound of 10 isolated spikes rose 4 to 18 ms after
  the time scheduled for their column (the measurement resolves one 17 ms frame); the
  whole-session view sounds throughout; over 123 scroll steps a column's cell values did not
  change. An independent review found that the first version went silent whenever the view
  stood still (the whole-session view among them), which the self-test then could not see; the
  scheduling now runs every frame and the transport checks above were added. Not checked by
  ear by the author of the code. Only Chrome on Windows was run; Safari and Firefox are untested.

## Data files

`data/` holds everything the viewer loads:
- `meta.json`: per-neuron table (cluster id, probe, shank, position, session rate), shank
  segments, colormaps, the 1 s population rate for the overview strip.
- `rastermap.json`: the Rastermap order (from `rastermap/`, via `prep_rastermap.py`).
- `manifest.json` and `spikes_ticks.gz.NN`, `spikes_rows.gz.NN`: all 60,997,315 spikes,
  packed by `pack_data.py`. Spike times are stored as differences from the previous spike
  (one byte each, since spikes are time-sorted) and neuron ids as two byte planes, each stream
  gzip-compressed and split into parts under 45 MB (GitHub refuses files over 100 MB). The
  browser decompresses them while downloading and checks the result against checksums of the
  original arrays; `pack_data.py` checks the round trip is exact.

To rebuild from the original recording (on the Windows desktop, where the inputs are):
`python prep_data.py` writes the raw `data/ticks.u32` and `data/rows.u16` (not in git), then
`python pack_data.py data` packs them. The viewer and `reference.py` use the raw files when
present, and the packed ones otherwise.

## Files

| File | Role |
|---|---|
| `index.html`, `viewer.js` | The viewer (no external libraries) |
| `sonify.js` | Plays the view as sound (Web Audio API); `?sonifytest=1` |
| `serve.py` | Local server for this folder (Windows or macOS) |
| `prep_data.py` | ALF + movie export to raw `ticks.u32`, `rows.u16`, and `meta.json` |
| `pack_data.py` | Packs the raw spike files for the web; `load_spikes()` reads either form |
| `prep_rastermap.py` | Copies the Rastermap ordering to `data/rastermap.json` (checks it is a permutation of the rows) |
| `prep_orderings.py` | Copies the further orderings to `data/orderings.json` with their menu labels (same checks) |
| `reference.py`, `make_refs.py` | Reference implementation and test-case writer |
| `rastermap/` | How the Rastermap order was computed and checked (its own README) |
| `bursts/` | The burst definition: detector, diagnostic figures, and literature notes (its own README) |
| `plan/` | The original design plan, example stills and the script that made them (paths there point to the Windows desktop) |
| `.nojekyll` | Tells GitHub Pages to serve the files as they are |

## Known limits and next steps

- GPU memory: about 610 MB (the spikes twice, by time and by neuron), plus about 370 MB of
  ordinary memory. Fine on the TITAN V and the M3 Max; a laptop with a small integrated GPU
  may fail to start or fall back to the slower counting method. Keeping only the search
  method would roughly halve the GPU memory.
- The online copy is public: anyone with the address can view and download the data.
- Each visit downloads 115 MB (the browser may cache it). GitHub Pages has a soft limit of
  100 GB of traffic a month, about 850 visits.
- Not yet tried on Windows: the Alt-key readout (releasing Alt alone can focus the browser's
  menu on Windows; the viewer blocks the default action but this is untested), the search
  counting method, and the three vertical layouts. Run `?test=1` there once.
- Future ideas noted: behavior event markers; sonification. (Rastermap order is now a mode.)
