"""Behavioral traces of FD_016 2026-11-13/1 for the raster viewer: the top three singular
components of the video's motion energy, and wheel velocity.

python make_fd016_traces.py <session folder> <motion_svd.npz> <timing_check.json> <viewer data folder>

Video. motion_svd.npz comes from extract_frames.py + motion_svd.py (crop chosen with
pick_crop.py). The video holds 15,969 frames, although the camera was triggered 169,919 times
(70 Hz, the whole 40 min session): it stops about 4 min into the session, so these traces
cover 3.5 to 232 s only. Frame times: frame 0 at the first trigger pulse, then one frame every
1 / rate, with rate = the fitted 69.93 frames/s from check_timing.py rather than the trigger's
70.00 (assuming "frame i at trigger i" puts the end of the video about 0.2 s early relative to
the firing rates). The times therefore rest on a fit to the neural data: good to about 50 to
75 ms, probably worst (traces up to about 50 ms early) around 115 to 170 s; see the notes in
check_timing.py.
Each component is centered on its median, which is its level while the mouse is still, and
divided by its standard deviation over the video: zero means rest, and the unit is that sd.

Wheel sign: positive when the counter counts up; which way the wheel then turns was not
established.

Wheel. rotaryEncoder in Timeline is the position counter (unsigned 32 bit, so it wraps at zero).
Velocity = difference between successive 2 kHz samples, smoothed with a Gaussian (sd 50 ms),
in encoder counts per second, sampled every 10 ms. The wheel barely turns in this session
(169 single-count steps in 40 min).
"""
import json
import sys

import numpy as np
from scipy import ndimage

from check_timing import TL_FS, trigger_times
from prep_traces import write_traces

WHEEL_SD = 0.05          # s
WHEEL_STEP = 20          # Timeline samples per output sample (10 ms)
# Not shown when the viewer opens (N. Steinmetz, 2026-11-25): the video covers under 4 of the
# 40 min and its frame times are a fit, and the wheel hardly moves in this session.
SHOW = False
WHEEL_YLIM = [-40.0, 40.0]   # counts/s; one single-count step peaks at 8, the one large movement (75 s) at 316 is clipped


def wheel_velocity(session):
    """(t0, dt, counts/s): smoothed velocity of the wheel from the Timeline position counter."""
    pos = np.load(f"{session}/rotaryEncoder.raw.npy").ravel()
    d = np.diff(pos)
    d = (d + 2.0 ** 31) % 2.0 ** 32 - 2.0 ** 31          # the counter is unsigned 32 bit
    assert np.abs(d).max() <= 1000, "implausible jump of the wheel counter"
    v = d * TL_FS                                        # counts/s; sample k is at (k + 0.5) / TL_FS
    # direct convolution: exactly zero wherever the wheel has been still for 4 sd (an FFT
    # convolution would leave round-off noise there)
    vs = ndimage.gaussian_filter1d(v, WHEEL_SD * TL_FS, mode="constant", truncate=4.0)
    first = WHEEL_STEP // 2 - 1
    out = vs[first::WHEEL_STEP]
    # the smoothed velocity integrates to the distance moved (edges aside)
    assert abs(out.sum() * WHEEL_STEP / TL_FS - d.sum()) < 0.5, (out.sum() * WHEEL_STEP / TL_FS, d.sum())
    return (first + 0.5) / TL_FS, WHEEL_STEP / TL_FS, out, int(np.count_nonzero(d)), float(d.sum())


def main():
    session, svd_file, timing_file, data_dir = sys.argv[1:5]
    z = np.load(svd_file)
    tim = json.load(open(timing_file))
    rate = tim["fitted_rate_hz"]
    trig = trigger_times(session)
    assert abs(trig[0] - tim["first_trigger_s"]) < 1e-9
    temporal = z["temporal"][:, :3].astype(np.float64)
    zs = (temporal - np.median(temporal, axis=0)) / temporal.std(axis=0)      # zero = the resting level
    t0 = trig[0] + 0.5 / rate                              # sample j: between frames j and j + 1
    crop = [int(a) for a in z["crop"]]
    print(f"video: {zs.shape[0]} motion samples, {t0:.4f} to {t0 + (zs.shape[0] - 1) / rate:.4f} s at {rate:.2f} per s; "
          f"under 'frame i at trigger i' the last sample would be at {(trig[zs.shape[0] - 1] + trig[zs.shape[0]]) / 2:.4f} s")
    note = (f"crop x {crop[0]} to {crop[2]}, y {crop[1]} to {crop[3]} of eye.mp4; frame times from the first camera "
            f"trigger at a fitted {rate:.2f} frames/s (see traces/check_timing.py), good to about 0.05 to 0.075 s; video ends at 232 s")
    traces = [dict(name=f"motion{k + 1}", label=f"Video motion {k + 1}", unit="sd", t0=t0, dt=1 / rate, values=zs[:, k], show=SHOW,
                   description=f"Singular component {k + 1} of the video's motion energy "
                               f"({100 * float(z['var_frac'][k]):.1f}% of its variance), minus its median (rest), "
                               f"in units of its standard deviation; " + note)
              for k in range(3)]
    wt0, wdt, wv, n_steps, net = wheel_velocity(session)
    print(f"wheel: {n_steps} counter steps in the session, net {net:.0f} counts; peak |velocity| {np.abs(wv).max():.1f} counts/s")
    traces.append(dict(name="wheel", label="Wheel velocity", unit="counts/s", t0=wt0, dt=wdt, values=wv, ylim=WHEEL_YLIM, show=SHOW,
                       description=f"Wheel velocity from the Timeline position counter: sample-to-sample difference, Gaussian "
                                   f"smoothing (sd {WHEEL_SD * 1e3:.0f} ms), encoder counts per second, positive when the counter "
                                   f"counts up; {n_steps} single-count steps in the whole session"))
    write_traces(data_dir, traces, note="FD_016 2026-11-13/1, made by traces/make_fd016_traces.py")
    np.savez(svd_file.replace(".npz", "_traces.npz"), t0=t0, dt=1 / rate, z=zs, wheel_t0=wt0, wheel_dt=wdt, wheel=wv)


if __name__ == "__main__":
    main()
