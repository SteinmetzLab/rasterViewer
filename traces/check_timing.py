"""Which Timeline time does each video frame belong to?

The camera is triggered by Timeline's 70 Hz clock output, recorded as cameraTrigger. The video
of FD_016 2026-11-13/1 has far fewer frames (15,969) than there are trigger pulses (169,919), so
the assumption "frame i was taken at trigger pulse i" needs a test against something with known
timing. The mouse's movements go with changes in firing across the brain (Stringer et al. 2019),
and the spike times are on the Timeline clock, so this script asks at which frame rate and start
time the video's motion energy lines up best with the population firing rates.

Model: frame i at time T + i / rate; the motion-energy sample between frames j and j + 1 at
T + (j + 0.5) / rate. For each rate the motion energy is put on a 20 ms grid and correlated with
the summed firing rate of each shank at every possible T in the session (sliding, normalized
correlation); the score is the root mean square over shanks of the correlation. A finer scan
near the trigger rate then fits the rate and a constant offset d, and the lag of firing behind
motion is estimated in 8 successive segments, which shows whether the timing drifts.

Result for FD_016 2026-11-13/1: over all rates from 20 to 100 Hz and all start times in the
session, the video lines up best at the trigger rate starting at the first trigger. Under
"frame i at trigger i" the lag grows from about 20 to about 240 ms through the video; at
69.93 frames/s it stays within 20 to 70 ms. make_fd016_traces.py uses the fitted rate.

What that does and does not establish (an independent review repeated the analysis with other
neural signals, filters and a synthetic control, 2026-11-25):
- the drift is real and not made by the analysis; the rate is 69.93 Hz with a rough 95%
  interval of 69.91 to 69.945 (block bootstrap), that is 15 frames (13 to 20) fewer than
  triggers over the video;
- it does not say why. Either frames were dropped, or the camera ran free at a 14.3 ms frame
  period (which is 69.93 Hz) from the first trigger on. Both give the same frame times;
- the drift is probably not uniform: lags by block look like steps (about 4 frames near 24 s,
  3 near 117 s, 2 near 71 and 145 s, singles elsewhere), which a staircase fits somewhat better
  out of sample than a line. That favors frames lost in bursts, but movement-dependent neural
  latency could look the same. With one rate, the traces may sit 40 to 50 ms early relative to
  the spikes around 115 to 170 s;
- a constant offset (frame 0 at trigger 1 to 3, or exposure delay) cannot be told from neural
  latency this way.
A frame counter or drop log from the acquisition software would settle it.

(A first attempt used the photodiode and the brightness of the video; the screen's light does
not show in this infrared video, so that had no signal.)

python check_timing.py <motion_svd.npz> <session folder> <viewer data folder> <out.json>
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

TL_FS = 2000.0
BIN = 0.02               # s: bins of the coarse scan
FINE_BIN = 0.005         # s: bins of the fine scan


def trigger_times(session):
    """Times (s, Timeline clock) of the rising edges of the camera trigger."""
    ct = np.load(f"{session}/cameraTrigger.raw.npy").ravel()
    thr = 0.5 * (np.percentile(ct, 1) + np.percentile(ct, 99))
    i = np.flatnonzero((ct[1:] > thr) & (ct[:-1] <= thr))
    # the edge lies between samples i and i + 1 (it is faster than a sample, so no finer estimate)
    return (i + 0.5) / TL_FS


def shank_rates(data_dir, n_bins, bin_s=BIN):
    """(n_groups, n_bins) spike counts per bin of each probe x shank, from the viewer's data."""
    sys.path.insert(0, str(Path(data_dir).resolve().parent))
    from pack_data import load_spikes
    meta = json.load(open(Path(data_dir) / "meta.json"))
    ticks, rows = load_spikes(Path(data_dir))[:2]
    grp = np.asarray(meta["rows"]["probe"]) * 4 + np.asarray(meta["rows"]["shank"])
    b = (ticks.astype(np.int64) // int(round(bin_s * meta["fs"]))).astype(np.int64)
    keep = b < n_bins
    n_grp = int(grp.max()) + 1
    flat = grp[rows[keep]].astype(np.int64) * n_bins + b[keep]
    return np.bincount(flat, minlength=n_grp * n_bins).reshape(n_grp, n_bins).astype(np.float64)


def highpass(x, cutoff=0.05, bin_s=BIN):
    sos = signal.butter(2, cutoff, "highpass", fs=1 / bin_s, output="sos")
    return signal.sosfiltfilt(sos, x, axis=-1)


def sliding_corr(templ, x):
    """Normalized correlation of templ (m,) with every length-m stretch of x (n,): (n - m + 1,)."""
    m = templ.size
    t = (templ - templ.mean()) / templ.std()
    dot = signal.fftconvolve(x, t[::-1], mode="valid")
    c1 = np.concatenate([[0], np.cumsum(x)])
    c2 = np.concatenate([[0], np.cumsum(x * x)])
    s1, s2 = c1[m:] - c1[:-m], c2[m:] - c2[:-m]
    sd = np.sqrt(np.maximum(s2 / m - (s1 / m) ** 2, 1e-12))
    return dot / (m * sd)                                # t has zero mean, so x's mean drops out


def main():
    svd_file, session, data_dir, out = sys.argv[1:5]
    z = np.load(svd_file)
    me = np.log(z["total_me"].astype(np.float64) + 1e-3)   # mean motion energy over the crop
    n_me = me.size
    trig = trigger_times(session)
    rate0 = (trig.size - 1) / (trig[-1] - trig[0])
    n_bins = int(np.floor(trig[-1] / BIN)) + 500
    R = highpass(shank_rates(data_dir, n_bins))
    print("spike counts per 20 ms:", R.shape, " trigger rate", rate0, "Hz, first trigger", trig[0], "s")

    def score_rate(rate):
        dur = n_me / rate
        tg = np.arange(BIN / 2, dur, BIN)                 # grid from the video's start
        m = np.interp(tg, (np.arange(n_me) + 0.5) / rate, me)
        m = highpass(m)
        cc = np.stack([sliding_corr(m, r) for r in R])    # (28, n_lags); lag k: video starts at k * BIN
        return np.sqrt((cc ** 2).mean(axis=0)), cc

    rates = np.unique(np.round(np.concatenate([np.arange(20, 100.01, 0.5), rate0 + np.arange(-1, 1.0001, 0.02)]), 4))
    best = []
    for rt in rates:
        s, _ = score_rate(rt)
        k = int(s.argmax())
        best.append((rt, s[k], k * BIN))
    best = np.array(best)
    o = np.argsort(-best[:, 1])
    print("best rates (Hz, score, start of video in s):")
    for i in o[:8]:
        print("   %.2f  %.4f  %.2f" % tuple(best[i]))
    s0, cc0 = score_rate(rate0)
    k0 = int(round(trig[0] / BIN))
    order = np.sort(s0)[::-1]
    print(f"at the trigger rate: best start {s0.argmax() * BIN:.2f} s (score {s0.max():.4f}); "
          f"score at the first trigger ({trig[0]:.2f} s): {s0[k0]:.4f}; next best elsewhere "
          f"(more than 2 s away): {np.delete(s0, np.arange(max(0, s0.argmax() - 100), s0.argmax() + 100)).max():.4f}; "
          f"median over starts {np.median(s0):.4f}")
    # finer, near the trigger rate: frame i at (first trigger) + d + i / rate. Spike counts in
    # 5 ms bins smoothed with a Gaussian (sd 20 ms); the motion energy interpolated at bin centers.
    FB = FINE_BIN
    nb = int((trig[0] + n_me / rate0 + 30) / FB)
    Rf = shank_rates(data_dir, nb, FB)
    gw = signal.windows.gaussian(81, 4)
    Rf = highpass(signal.fftconvolve(Rf, (gw / gw.sum())[None], mode="same"), 0.1, FB)
    tb = (np.arange(nb) + 0.5) * FB

    def score(t_me, t_a=None, t_b=None, lag=0.0):
        """Root mean square over shanks of the correlation of motion energy at times t_me with
        the firing rates lag seconds later, between t_a and t_b (2 s inside the video by default)."""
        i0 = int(((t_me[0] + 2) if t_a is None else t_a) / FB)
        i1 = int(((t_me[-1] - 2) if t_b is None else t_b) / FB)
        m = np.interp(tb[i0:i1], t_me, me)
        m = m - m.mean()
        k = int(round(lag / FB))
        seg = Rf[:, i0 + k:i1 + k]
        c = ((seg - seg.mean(1, keepdims=True)) * m).sum(1) / (seg.std(1) * m.std() * m.size)
        return float(np.sqrt((c ** 2).mean()))

    frates = np.arange(rate0 - 0.20, rate0 + 0.0201, 0.01)
    ds = np.arange(-0.10, 0.3001, 0.005)
    S = np.array([[score(trig[0] + d + (np.arange(n_me) + 0.5) / r) for d in ds] for r in frates])
    ia, ib = np.unravel_index(S.argmax(), S.shape)
    j0, k0d = int(np.argmin(np.abs(frates - rate0))), int(np.argmin(np.abs(ds)))
    print(f"fine scan: best rate {frates[ia]:.2f} Hz with d = {ds[ib] * 1e3:.0f} ms (score {S[ia, ib]:.4f}); "
          f"at the trigger rate: best d {ds[S[j0].argmax()] * 1e3:.0f} ms (score {S[j0].max():.4f}), d = 0: {S[j0, k0d]:.4f}")

    # does the lag of firing behind motion drift through the video? 8 segments, two models
    lags = np.arange(-0.5, 0.5001, FB)
    seg_lags = {}
    for label, t_me in (("trigger", (trig[:n_me] + trig[1:n_me + 1]) / 2),
                        ("fitted", trig[0] + (np.arange(n_me) + 0.5) / frates[ia])):
        edges = np.linspace(t_me[0] + 2, t_me[-1] - 2, 9)
        out_l = []
        for x, y in zip(edges[:-1], edges[1:]):
            sc = [score(t_me, x, y, l) for l in lags]
            out_l.append((float((x + y) / 2), float(lags[int(np.argmax(sc))]), float(max(sc))))
        seg_lags[label] = out_l
        print(f"lag of firing behind motion, model '{label}', by segment (ms):", [int(round(l * 1e3)) for _, l, _ in out_l])
    res = dict(n_motion_samples=int(n_me), n_triggers=int(trig.size), trigger_rate_hz=float(rate0),
               first_trigger_s=float(trig[0]), best_rate_hz=float(best[o[0], 0]), best_rate_score=float(best[o[0], 1]),
               best_rate_start_s=float(best[o[0], 2]), start_at_trigger_rate_s=float(s0.argmax() * BIN),
               score_at_trigger_rate=float(s0.max()), score_at_first_trigger=float(s0[k0]),
               median_score_over_starts=float(np.median(s0)),
               fitted_rate_hz=float(frates[ia]), fitted_d_s=float(ds[ib]), fitted_score=float(S[ia, ib]),
               score_trigger_rate_d0=float(S[j0, k0d]), segment_lags=seg_lags,
               frames_missing_estimate=float(n_me * (rate0 / frates[ia] - 1)))
    json.dump(res, open(out, "w"), indent=1)
    np.savez(out.replace(".json", ".npz"), rates=best[:, 0], rate_score=best[:, 1], rate_start=best[:, 2],
             start_score=s0, bin=BIN, cc_at_rate0=cc0.astype(np.float32), fine_rates=frates, fine_d=ds, fine_score=S)


if __name__ == "__main__":
    main()
