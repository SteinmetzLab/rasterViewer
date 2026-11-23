# Raw-voltage underlay (prototype, 2026-11-23)

With the neurons plotted by depth and a window of 0.5 s or less, the viewer draws the recorded
voltage of every site under the spikes. This is a local prototype on one stretch of one probe:
62 s of probe 1 (`p1_imec0`, 4 shanks, 1,536 sites) around t = 1100 to 1160 s. The data are
not in the repository and not on the public site.

## Running it

    python serve.py --raw Y:/temp/fd016_raw_underlay

then open `http://localhost:8765/#ymode=depth&t0=1130&t1=1130.05&n0=0&n1=2880`. The yellow mark
under the overview strip shows where raw data exist. On a Mac, give `--raw` the path at which
the lab server's `temp/fd016_raw_underlay` folder is mounted.

Keys (only when neurons are plotted by depth):

| Key | Does |
|---|---|
| V | next stored version: lossless, lossy WavPack, lossy with 3 bits dropped |
| D | difference between the current lossy version and lossless |
| [ and ] | halve / double the gray range (starts at ±100 µV) |
| O | spike marks at 100%, 60%, 0% opacity |
| U | hide / show the voltage |

## What is drawn

- **Processing** (in a worker, `raw_worker.js`; the same in every version): per site, a 300 Hz
  high-pass (Butterworth, order 3, run forward and backward); then, per sample, the median over
  the shank's working sites is subtracted. Dead sites (from the sorter's `provenance.json`) are
  flat. This is simpler than what the sorter saw (no correction for the sites' sampling
  delays, one median over the shank and not one per digitizer group).
- **Pixels**: every screen pixel shows the sample of largest magnitude among the samples and
  sites it covers. Black is negative.
- **Depth**: a site at height y occupies the same 15 µm band as a neuron at that height. The two
  sites of a row (x = 27 and 59 µm) share the band: left site in the upper half, right site in
  the lower half. So the two columns of sites appear interleaved, as in a trace viewer.
- **Time**: spike times in the viewer are on the behavioral rig's clock, where they were put by
  piecewise-linear interpolation between sync edges. That map wobbles around a straight line by
  up to 0.3 ms (sd 0.12 ms), because the rig samples the sync signal at 2 kHz. The underlay
  sends viewer time through the same map (recovered from the spikes: every sorted spike has a
  sample number and a time), so a spike mark lands on its waveform, to within about one
  sample (spike times are rounded to whole ticks and drawn as cells). It also means that spike times of different probes carry that wobble relative to each
  other in the data set itself.
- The spike marks are the viewer's usual cells, drawn over the voltage. Cells without spikes are
  transparent there.

## Stored versions

All three are Zarr (version 2) arrays, one per shank, 384 sites x samples, in chunks of 0.2 s,
with a delta filter along time and zlib. A browser decodes that with nothing but its built-in
`DecompressionStream`.

| Version | What it holds | Smaller than the recorded int16 by | Error, sd | Largest error |
|---|---|---|---|---|
| `lossless` | the recorded samples | 2.67× (WavPack lossless would give 3.46×) | none | none |
| `wv2.25` | samples after lossy WavPack at 2.25 bits per sample, decoded | 7.11× (as WavPack) | 5.8 µV | 73 µV |
| `bits3` | samples with the 3 lowest bits dropped (steps of about 24 µV) | 5.97× | 7.0 µV | 10.6 µV |

(All four shanks, 63 s; errors over all samples before any filtering. `manifest.json` has the
numbers.) For scale, the noise of a site after the display processing is about 14 µV
(shank 0), so the lossy versions raise it to about 15.4 µV (WavPack, 7% more) and 16.0 µV
(bits dropped, 11% more). On shank 0 (10 s, after the high-pass):

- WavPack's error grows a little with the signal: sd 5.5 µV where the voltage is under 30 µV,
  6.9 µV where it is over 100 µV. The error of dropped bits is the same everywhere (6.8 µV).
- Neither shifts spike troughs: at 58,385 troughs below -80 µV the mean error is -0.6 µV
  (WavPack) and 0.0 µV (bits dropped), sd 6.3 and 6.8 µV.
- Not measured: the effect on spike sorting. Buccino et al. report a slight loss of low-amplitude
  units on NP2.0 at this WavPack level.

The WavPack level is the strongest one tested by Buccino et al. (eLife 2026,
doi 10.7554/eLife.110170). No WavPack decoder for browsers was found, so `wv2.25` stores the
decoded samples (which is what such a decoder would return) and is no smaller than lossless on
disk; the manifest records what WavPack itself needs. `bits3` is the real format: what is on
disk is what a browser would download.

## Files

| File | What |
|---|---|
| `make_raw_proto.py` | reads the lab's `.cbin` files, writes the arrays, the clock pairs, the manifest and the test references |
| `cbin.py` | reads pieces of an mtscomp `.cbin` |
| `validate_alignment.py` | checks from the raw data that each unit's spikes sit on its waveform, in time and in depth |
| `../raw_underlay.js`, `../raw_worker.js` | the viewer side |

## Checks

- `make_raw_proto.py` reads every array back and compares it with what it wrote.
- `http://localhost:8765/?rawtest=1` compares the browser's decoding (sum of a chunk of every
  array, 12 of 12 equal) and its processing (one tile of every version) with the Python
  reference (`process_tile`): identical to the last bit of the float32 values.
- The viewer's own `?test=1` is unchanged by the new code (120 runs, 0 mismatched cells).
- Reading of the `.cbin` files: every site with zero variance in the decoded data is on the
  sorter's list of dead sites (shanks 0 to 2: the lists are equal; shank 3: the sorter lists one
  more site, 153, which is not flat in this stretch).
- `validate_alignment.py`, all four shanks, 460 units with a spike-triggered average of 30 µV
  or more (results in `validate_alignment_sh<k>.json` next to the data):
  - Time: the trough of the average is 0 or 1 sample before the mapped spike time for 82% of
    units and within 2 samples for 90% (by shank: 94, 88, 88, 92%). With a straight-line clock
    in place of the pairs, the averages shrink to 54% of their amplitude (median), i.e. the
    sync wobble is real and the pairs remove it.
  - Depth: the site with the largest trough is within 15 µm of the unit's depth in the viewer
    for 84% of units (median distance about 2 µm; by shank: 86, 86, 79, 88%). For 37 units (8%)
    the largest trough within 150 µm is 60 µm or more away. On shank 0, where these were looked
    at, the units' own template peaks agree with their viewer depth, and the larger trough
    belongs to a bigger waveform that occurs at the same time; see "Known limits".

## Known limits

- **Motion.** The sorter ran on motion-corrected data, so unit depths are in corrected
  coordinates and the voltage here is not. On probe 1 the estimated motion is small (at most
  about 7 µm on shanks 0 to 2; 25 µm at its largest on shank 3), but other probes moved by up
  to 70 µm. A full version should shift the voltage by `sorting/<shank>/motion.npz`.
- **Several units on one waveform.** Some large spikes extend over 100 µm or more along the
  shank and carry the marks of two or more sorted units at different depths. Whether those are
  one neuron sorted more than once has not been examined.
- **One probe, 63 s, one stretch.** Other probes need their own clock pairs and arrays.
- **Tile edges.** A pixel column that straddles two tiles (one every 0.2 s) shows only the
  samples of the tile holding its central sample, and the very last sample of the data is not
  drawn. An independent check of the drawing (12 million pixels in 7 views, against the rule
  computed separately) found 0.05% of pixels different, at these seams and where a pixel edge
  falls within 0.0005 sample of a sample boundary.
- **Speed** (this desktop, data on the lab server): one shank and 50 ms in 0.45 s, the whole
  probe and 0.3 s in 1.3 s, 0.5 s in 2.0 s; instant once loaded.
- **WavPack in a browser** needs a decoder compiled to WebAssembly; none was found ready-made.
