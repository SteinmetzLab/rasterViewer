# Behavioral traces

The viewer draws regularly sampled series under the raster when `data/traces.json` exists
(`traces.js`; see "Behavioral traces" in the main README). This folder holds the writer for
those files and the scripts that made the traces of FD_016 2026-11-13/1.

## Writing traces for a data set

```python
from prep_traces import write_traces
write_traces("../data", [
    dict(name="wheel", label="Wheel velocity", unit="cm/s", t0=0.005, dt=0.01, values=v),
])
```

Each trace: `name` (letters, digits, `_`; used in the URL), `label` and `unit` (shown beside the
lane), `t0` and `dt` (sample i belongs to time `t0 + i * dt` seconds on the viewer's clock, the
clock of the spike times), `values` (NaN where data are missing), and optionally `ylim`, the
fixed vertical range of the lane. Without `ylim` the range is the 0.5th to 99.5th percentile,
widened to include zero. The viewer shades between zero and the trace, so give it series for
which zero means something (rest, baseline, the mean).

`python prep_traces.py ../data` checks the files of an existing list.
`python test_prep_traces.py` tests the writer.

## The traces of FD_016 2026-11-13/1

Four traces: the top three singular components of the motion energy of the face video, and
wheel velocity.

| Step | Script | Output |
|---|---|---|
| Choose the crop | `pick_crop.py` (a window: drag a rectangle, click OK) | `crop.json` |
| Read the video once | `extract_frames.py` | cropped frames averaged in 4 x 4 pixel blocks |
| Motion energy and its SVD | `motion_svd.py <frames> <out> 10 30` | components, masks, variance |
| Frame timing | `check_timing.py` | fitted frame rate, lag by segment |
| Trace files | `make_fd016_traces.py` | `data/traces.json`, `data/trace_*.f32` |
| Figure | `fig_motion.py` | crop, masks, time courses, variance |

**Motion energy** (after Stringer et al. 2019, Science 364:eaav7893, as in Facemap): the
absolute difference between consecutive frames per binned pixel; the mean image is subtracted
and the singular value decomposition of the (time x pixel) matrix gives the components. The
crop is x 78 to 698, y 242 to 702 of the 1280 x 1024 frame (face, ear, forepaw; chosen by
N. Steinmetz). Components 1 to 3 hold 21.7%, 4.5% and 2.8% of the variance: overall motion
(strongest at the whisker pad; correlation 0.93 with the total motion energy), whisker pad
against body and ear, and mostly ear and forepaw. Signs are set so that each mask sums to a
positive number. For the viewer each component is centered on its median (its level while the
mouse is still) and divided by its standard deviation.

**Keyframe comb.** The video is H.264 with a full image every 30 frames. The frame pair at
each full image has 31% more motion energy than the others when nothing moves, which would
show as a 2.3 Hz rhythm in the traces. `motion_svd.py` subtracts, per pixel and per phase of
the 30 frames, the median excess; 3.5% remains at that phase.

**What the video covers, and when its frames were taken.** The camera was triggered at 70 Hz
for the whole 40 min session (169,919 pulses on Timeline), but `eye.mp4` holds 15,969 frames
and was last written about 5 min into the session. The reason is not known. The traces
therefore cover 3.5 to 232 s only.

Since frames and triggers cannot be paired by count, `check_timing.py` compares the motion
energy with the firing rates of the 28 shanks (spike times are on the Timeline clock):

- Over all frame rates from 20 to 100 Hz and all start times in the session, the best match is
  the trigger rate starting at the first trigger pulse (score 0.26; at most 0.09 at any other
  start time).
- With frame i at trigger i, firing lags motion by about 20 ms in the first eighth of the
  video and 240 ms in the last. At 69.93 frames/s the lag stays within 20 to 70 ms. The trace
  files use 69.93 frames/s from the first trigger.
- An independent review (other neural signals and filters, a synthetic control with known
  truth, a block bootstrap) confirmed the drift and gave a rough 95% interval of 69.91 to
  69.945 frames/s, that is 13 to 20 frames fewer than triggers over the video.
- Why is open: frames dropped by the acquisition, or a camera running free at a 14.3 ms frame
  period (69.93 Hz) rather than on the trigger. The lag by block looks more like steps than a
  line, which favors frames lost in bursts; in that case the single fitted rate puts the
  traces up to about 50 ms early around 115 to 170 s.
- So the video traces are placed to about 50 to 75 ms, by a fit to the neural data. Do not
  use them to measure latencies between movement and firing. A frame counter or drop log from
  the camera software would settle the timing.

**Wheel.** `rotaryEncoder` in Timeline is the position counter (unsigned 32 bit, it wraps at
zero). Velocity is the difference between successive 2 kHz samples, smoothed with a Gaussian
(sd 50 ms), in encoder counts per second, stored every 10 ms; positive when the counter counts
up (the direction of the wheel for that was not established, and counts were not converted to
distance). The wheel barely turns in this session: 169 single-count steps in 40 min, net 59
counts, the only sizable movement at 75 s (peak 316 counts/s; the lane is fixed at +-40, so it
is clipped there).

## Checks

- `?tracetest=1` in the viewer: 25 of 25 checks (files against recorded counts and sums; exact
  column summaries; drawing of known series; placement).
- `motion_svd.py` asserts that the factors are orthonormal and satisfy the SVD equations.
- `make_fd016_traces.py` asserts that the smoothed wheel velocity integrates to the distance
  moved; `prep_traces.py` reads every file back.
- Review by a separate agent (2026-11-25): pyramid summaries against a plain loop on 1.5
  million ranges, the wheel unwrapping and time stamps, the trigger detection, the comb
  removal, and the timing analysis above. Its findings were applied (readout column, redraw
  while Sonify plays, one-color shading for lanes that do not cross zero, exact zeros in the
  wheel trace, a least lane height, a more robust default range, zero at rest for the motion
  components).

## Limits

- Values beyond a lane's range are clipped without a mark (the readout shows them): 0.5% of
  the samples of component 1.
- A sample lying exactly on a column edge may fall on either side by rounding; no sample is
  skipped or counted twice.
- One recording; the crop and the binning (4 x 4) were not varied.
