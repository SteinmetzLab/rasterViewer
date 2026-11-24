# Raw-voltage underlay

With the neurons plotted by depth, a window of 0.5 s or less and a few shanks in view, the
viewer draws the recorded voltage of every site under the spikes: all 7 probes (28 shanks,
10,752 sites), the whole synced period (1.97 to 2,430.81 s). The data are hosted in a
Cloudflare R2 bucket, which is public: anyone with the address can download them.

This folder holds the code that made the hosted data and the checks on it. The viewer side is
`../raw_underlay.js` and `../raw_worker.js`.

## What is drawn

- **Processing** (in a worker, `raw_worker.js`): per site, a 300 Hz high-pass (Butterworth,
  order 3, run forward and backward); then, per sample, the median over the shank's working
  sites is subtracted. Dead sites (from the sorter's `provenance.json`) are flat. This is
  simpler than what the sorter saw (no correction for the sites' sampling delays, one median
  over the shank and not one per digitizer group, no interpolation of bad sites).
- **Pixels**: every screen pixel shows the sample of largest magnitude among the samples and
  sites it covers. Black is negative.
- **Depth**: a site at height y occupies the same 15 µm band as a neuron at that height. The two
  sites of a row (x = 27 and 59 µm) share the band: left site in the upper half, right site in
  the lower half. So the two columns of sites appear interleaved, as in a trace viewer.
- **Time**: spike times in the viewer are on the behavioral rig's clock, where they were put by
  piecewise-linear interpolation between sync edges. That map wobbles around a straight line by
  up to about 0.5 ms, because the rig samples the sync signal at 2 kHz. The underlay sends
  viewer time through the same map, one per probe (`clock_knots` in `make_raw.py`), so a spike
  mark lands on its waveform to within about one sample.
- **Motion**: the sort ran on data that SpikeInterface had motion-corrected with a DREDGE
  estimate (7 depth windows, 1 s bins; `sorting/<shank>/motion.npz`), so the depth of every
  unit is in corrected coordinates: corrected(y) = recorded(y + d). The stored voltage is as
  recorded, and the viewer shifts it when drawing: a pixel at depth y shows the sites at
  y + d, with d the estimate for that second at that depth (linear between the depth
  windows, constant beyond the outer ones). The estimate is under 5 µm in 98.8% of all
  (shank, second, depth window) of the synced period and over 15 µm in 0.05%; it reaches 52 µm
  on probe 2, shank 0 and 59 µm on probe 5, shank 3, in short excursions. Where the shift
  moves the end of the shank into view, nothing is drawn (mid gray). Key M switches the shift
  off.
- The spike marks are the viewer's usual cells, drawn over the voltage. Cells without spikes are
  transparent there.

A view is eligible when (shanks in view) × (0.2 s tiles the window can touch) is at most 24:
12 shanks at 0.2 s or less, 8 at 0.4 s, 6 at 0.5 s. Otherwise the status line says what to
zoom.

## The hosted data

Written by `make_raw.py`; layout under `https://pub-6f7885d91fac46b4a0aad6e22005454a.r2.dev/v1/`:

| File | Content |
|---|---|
| `manifest.json` | probes, shanks, dead sites, scaling, tile ranges, the filter's coefficients, test references |
| `p<probe>/sh<shank>/0.<tile>` | one tile: 0.2 s (6,000 samples) of the shank's 384 sites |
| `p<probe>/sh<shank>/.zarray` | makes each shank's folder a Zarr (version 2) array, 384 × samples |
| `clock_p<probe>.f64` | knots (probe sample, viewer time) of the probe's clock map, float64 pairs |
| `motion_p<probe>.f32` | the motion estimate, float32 (shank, 1 s bin, depth window) |
| `test/` | reference output of the display processing for three tiles |

Probes are numbered 0 to 6 in the order `p1_imec0` to `p1_imec3`, `p2_imec0` to `p2_imec2`
(the viewer shows them as Probe 1 to 7). Tile T holds samples T × 6,000 to T × 6,000 + 5,999 of
the probe's recording. Only the tiles of the synced period are stored (793 to about 12,940);
a Zarr reader returns zeros outside.

- **Samples**: the recorded int16 samples with their 2 lowest bits dropped
  (stored = sample >> 2; the viewer uses stored × 4 + 1.5). NP2.0 Quad-Base probes digitize
  at 12 bits, 3.027 µV per count, so the stored step is 12.1 µV and the error is uniform
  within ±4.5 µV (sd 3.4 µV), which raises a site's noise of about 14 µV by about 3%.
  Nothing else is done to the samples: no filter, no reference, no motion correction.
- **Chunk format**: delta along time, then zlib (level 4). A browser decodes a chunk with
  its built-in `DecompressionStream("deflate")` and a running sum; zarr-python reads the
  arrays as they are.
- **Size**: SIZE_PLACEHOLDER
- **Cost and limits**: R2 charges $0.015 per GB per month for storage and nothing for
  downloads. The `r2.dev` address is meant by Cloudflare for development: it is rate-limited
  (hundreds of requests per second) and not cached at their edge. Tiles are served with a
  one-year browser cache, so a tile must never change under `v1/`; a changed format goes
  under a new prefix. The bucket only allows pages from `www.steinmetzlab.net`,
  `steinmetzlab.net` and `localhost:8765` / `localhost:8771` to read it from a browser.

## Running the viewer on other copies of the data

- `python serve.py --raw Y:/temp/fd016_raw_v1` serves a local copy as `rawdata/`, and a viewer
  on localhost uses it in place of the bucket.
- `?rawbase=<address>` reads from another address; `?rawbase=off` switches the underlay off.

## Making the data

    python make_raw.py meta    --out Y:/temp/fd016_raw_v1
    python make_raw.py convert --out Y:/temp/fd016_raw_v1 --keys KEYS.json
    python make_raw.py refs    --out Y:/temp/fd016_raw_v1
    python make_raw.py sync    --out Y:/temp/fd016_raw_v1 --keys KEYS.json

`convert` reads the lab's per-shank `.cbin` files one second at a time, writes the tiles to
the output folder and uploads each one; it can be stopped and restarted. `sync` compares the
bucket with what was written (name, size, MD5), uploads what is missing, and fails while any
tile of the manifest's range is absent. KEYS.json holds the bucket's credentials and is not in
this repository. The whole recording took TIME_PLACEHOLDER on one desktop (13 processes).

## Checks

- **Stored samples** (`verify_raw.py`): tiles read through zarr-python are compared with the
  same samples read from the `.cbin` files by the `mtscomp` package, neither being the
  converter's own code. VERIFY_PLACEHOLDER
- **Upload** (`make_raw.py sync`): SYNC_PLACEHOLDER
- **Clock** (`make_raw.py meta`, `clock_check.json`): the knots reproduce the time of every
  one of the 60,997,315 sorted spikes from its sample number to within 5 × 10⁻⁷ sample, on
  all 7 probes, and the map is continuous at every sync edge.
- **Browser against Python** (`?rawtest=1`, run against a local copy and against the
  bucket): clock and motion files of all 7 probes load; one chunk of each of the 28 shanks
  decodes to the stored sum; the display processing of three tiles equals the Python
  reference (`process_tile`) to the last bit of the float32 values; and the drawing of four
  views (one shank; a shank that had moved by 52 µm; four shanks across a tile edge and a
  change of the motion estimate; the shift switched off) is recomputed pixel by pixel in
  double precision: 0.001% to 0.06% of 1.9 million pixels differ, where a pixel edge falls on
  a sample or site boundary. The motion values and clock lookups the page uses were also
  compared with the sorter's files directly.
- **Spikes on their waveforms** (`validate_raw.py`, spike-triggered averages of the stored
  voltage at the viewer's spike times, units of 30 µV or more):

  | Shank, stretch | Units | Trough within 2 samples of the spike mark | Average with a straight-line clock (fraction of amplitude kept) |
  |---|---|---|---|
  | Probe 2, shank 0, 4 to 72 s | 49 | 94% | 0.54 |
  | Probe 2, shank 0, 238 to 317 s | 55 | 95% | 0.55 |
  | Probe 2, shank 0, 1,097 to 1,163 s | 52 | 94% | 0.56 |
  | Probe 5, shank 3, 1,558 to 1,642 s | 38 | 92% | 0.66 |

  (The prototype had found the same on probe 1: 90% of 460 units within 2 samples.)
- **Motion shift**:
  - The rule was compared with SpikeInterface's own code (`interpolate_motion_on_traces`,
    version 0.104.8; the sort ran 0.104.3): the same displacement to 10⁻¹⁵ µm at 20,000 points,
    and a test bump placed at y + d in the recorded data comes out at y after correction, for
    both signs of d.
  - On probe 2, shank 0, at 238 to 317 s, 13 units sat where the estimate averaged 5 µm or
    more (up to 31 µm). The center of their waveform, in the viewer's coordinates, is off from
    the unit's depth by 3.1 µm (median) with the shift, 12.6 µm without it, and 26.2 µm with
    the opposite sign. At 4 to 72 s (9 units): 3.9, 10.9 and 19.3 µm.
  - The same units compared between 238 to 317 s and 1,097 to 1,163 s, when the shank was
    near rest (`compare_stretches.py`): all 9 units whose estimate differed by 10 µm or more
    had their waveform move in the direction of the estimate, by 65% of it on average. Four
    moved by the estimated amount to within 4 µm; the five that moved less include the three
    units within 60 µm of the shank's tip, where the waveform runs off the end of the shank
    and the estimate is an extrapolation. So the estimate's sign is certain and its size is about
    right; the shift reproduces what the sorter was given, which is what puts a unit's
    waveform at the unit's depth.
  - Probe 5, shank 3, has a short excursion of up to 59 µm near 1,590 s. Only 4 units sit
    where it is large; for them the offset is 16 µm with the shift, 29 µm without, 42 µm with
    the opposite sign. The alignment there is better with the shift but not good.

## Known limits

- **Several units on one waveform.** Some large spikes extend over 100 µm or more along the
  shank and carry the marks of two or more sorted units at different depths, probably spikes
  of different neurons that coincided and that the sorter separated.
- **Unmarked spikes.** Only units that passed quality control are in the viewer: 49% of the
  spikes Kilosort assigned. Clear waveforms without a mark are usually spikes of a cluster
  that failed.
- **Tile edges.** A pixel column that straddles two tiles (one every 0.2 s) shows only the
  samples of the tile holding its central sample.
- **Motion between seconds.** The shift changes once per second of probe time, as it did for
  the sorter, so the voltage can step by a few micrometers there.
- **Speed.** A tile is about 1.2 MB to download and takes SPEED_PLACEHOLDER to process; each
  tile also needs the edges of its two neighbors in time, which are downloaded whole. Viewing
  0.5 s of one probe (4 shanks) for the first time downloads 25 to 30 MB.
- **Not measured**: the effect of dropping 2 bits on spike sorting. It was chosen for display.

## Files

| File | What |
|---|---|
| `make_raw.py` | the converter: clocks, motion, manifest, tiles, upload, test references |
| `cbin.py` | reads pieces of an mtscomp `.cbin` |
| `verify_raw.py` | stored samples against the `.cbin` files, through zarr-python and mtscomp |
| `validate_raw.py`, `compare_stretches.py` | spikes on their waveforms in time and depth; the motion shift |
| `make_raw_proto.py`, `validate_alignment.py` | the prototype of 2026-11-23 (63 s of one probe in three stored versions); `make_raw.py` uses its clock pairing and its reference of the display processing |

## How the stored version was chosen (prototype, 2026-11-23)

Three versions of 63 s of probe 1 were compared in the viewer (all four shanks; errors over
all samples before any filtering):

| Version | Smaller than the recorded int16 by | Error, sd | Largest error |
|---|---|---|---|
| lossless (delta + zlib) | 2.67× (WavPack lossless would give 3.46×) | none | none |
| lossy WavPack, 2.25 bits per sample | 7.11× | 5.8 µV | 73 µV |
| 3 lowest bits dropped (24 µV steps) | 5.97× | 7.0 µV | 10.6 µV |

WavPack's error has structure (it grows with the signal), the error of dropped bits is flat,
and no WavPack decoder for browsers was found. Dropping bits on these 12-bit probes: 1 bit
(6.1 µV steps) 3.33×, +0.6% noise; 2 bits (12.1 µV) 4.32×, +2.7%; 3 bits (24.2 µV) 5.97×,
+11% (measured on 80 tiles of probe 1). Two bits were chosen. The WavPack level is the
strongest one tested by Buccino et al. (eLife 2026, doi 10.7554/eLife.110170).
