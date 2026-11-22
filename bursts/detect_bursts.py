"""Burst detection on the viewer's spikes (per unit, from spike times in 30 kHz ticks).

A burst under criterion (silence S, max ISI M, min spikes K):
  a run of K or more consecutive spikes of one unit in which every ISI is <= M, and whose first
  spike is preceded by at least S of silence (no spike of that unit). S = 0 drops the silence
  requirement. Every spike in a qualifying run is a "burst spike"; its first spike is the burst
  onset. All other spikes are "tonic".

The thalamic low-threshold (LTS) burst convention comes from Lu, Guido & Sherman 1992
(J Neurophysiol 68:2185): ISIs <= 4 ms, first spike after >= 100 ms of silence, >= 2 spikes
(as stated in e.g. Reinagel & Reid 2000; Whitmire et al. 2016 Cell Rep). The silence matters:
T-type Ca2+ channels need sustained hyperpolarization to de-inactivate, and without it fast
tonic firing (e.g. cortical fast-spiking cells) also passes a short-ISI test. The first spike
of each unit's recording is treated as preceded by silence.

    python detect_bursts.py      # runs the self-tests below
"""
import numpy as np

FS = 30000

CRITERIA = {
    # name: (silence_ms, max_isi_ms, min_spikes)
    "LTS (Lu 1992): >=100 ms silence, ISI <=4 ms, >=2 spikes": (100.0, 4.0, 2),
    "relaxed LTS: >=50 ms silence, ISI <=6 ms, >=2 spikes": (50.0, 6.0, 2),
    "no silence: ISI <=8 ms, >=3 spikes": (0.0, 8.0, 3),
}


def by_unit(ticks, rows, n_units):
    """Spike order sorted by (unit, time), and each unit's start offset. Stable, so each unit
    stays in time order (ticks are time-sorted)."""
    order = np.argsort(rows, kind="stable")
    off = np.zeros(n_units + 1, np.int64)
    np.cumsum(np.bincount(rows, minlength=n_units), out=off[1:])
    return order, off


def detect(t_unit_sorted, rows_unit_sorted, silence_ms, max_isi_ms, min_spikes, fs=FS):
    """Classify spikes given arrays sorted by (unit, time). Returns (is_burst, is_onset, burst_id)
    as boolean / int arrays in the same order; burst_id is -1 for tonic spikes."""
    t = t_unit_sorted.astype(np.int64); r = rows_unit_sorted.astype(np.int64)
    n = t.size
    pre = np.full(n, np.iinfo(np.int64).max)            # ticks since the unit's previous spike
    same = r[1:] == r[:-1]
    pre[1:][same] = (t[1:] - t[:-1])[same]
    M = max_isi_ms * fs / 1000.0
    S = silence_ms * fs / 1000.0
    link = pre <= M                                      # joined to the previous spike
    run = np.cumsum(~link)                               # run index; a run starts where link is False
    start = np.flatnonzero(~link)
    length = np.diff(np.append(start, n))
    ok = (length >= min_spikes) & (pre[start] >= S)
    good = ok[run - 1]
    burst_id = np.where(good, run - 1, -1)
    onset = np.zeros(n, bool); onset[start[ok]] = True
    return good, onset, burst_id


def _selftest():
    ms = FS // 1000
    # unit 0: tonic spike, then an LTS burst after 200 ms silence (3 spikes, 3 ms apart),
    # then a fast train with no preceding silence (5 spikes 5 ms apart), then a lone doublet
    # 2 ms apart after 60 ms silence
    t0 = [0, 200 * ms, 203 * ms, 206 * ms, 230 * ms, 235 * ms, 240 * ms, 245 * ms, 250 * ms,
          310 * ms, 312 * ms]
    # unit 1: two spikes 1 ms apart at the very start (first spike counts as after silence)
    t1 = [5 * ms, 6 * ms]
    t = np.array(t0 + t1); r = np.array([0] * len(t0) + [1] * len(t1))
    b, o, _ = detect(t, r, 100, 4, 2)
    assert b.tolist() == [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1], b.astype(int).tolist()
    assert o.tolist() == [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    b, o, _ = detect(t, r, 50, 6, 2)      # the 60 ms-silence doublet now counts; the train's
    assert b.astype(int).tolist() == [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]  # ISI 5 ms but no silence
    b, o, _ = detect(t, r, 0, 8, 3)       # no silence: the LTS burst and the train
    assert b.astype(int).tolist() == [0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    # the LTS burst (206 ms) and the train (230 ms) are 24 ms apart, so they are separate runs
    # under ISI <= 8 ms, but both pass; check the ids differ
    _, _, ids = detect(t, r, 0, 8, 3)
    assert ids[3] != ids[4]
    # by_unit reorders interleaved spikes and keeps time order within units
    tt = np.array([1, 2, 3, 4, 5]); rr = np.array([1, 0, 1, 0, 1])
    order, off = by_unit(tt, rr, 2)
    assert tt[order].tolist() == [2, 4, 1, 3, 5] and off.tolist() == [0, 2, 5]
    print("detect_bursts self-tests passed")


if __name__ == "__main__":
    _selftest()
