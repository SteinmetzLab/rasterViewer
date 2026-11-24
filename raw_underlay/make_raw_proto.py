"""Prototype data for the raw-voltage underlay of the raster viewer: a stretch of one probe.

    python make_raw_proto.py --out Y:/temp/fd016_raw_underlay [--probe 0] [--t0 1100] [--dur 60]

Reads the lab's compressed raw files (one .cbin per shank, see cbin.py) and writes, per shank,
one Zarr (version 2) array per "variant", channels x samples, int16:

  lossless   the recorded samples
  wv2.25     the samples after lossy WavPack compression at 2.25 bits per sample (the strongest
             level tested by Buccino et al., eLife 2026, doi 10.7554/eLife.110170) and
             decompression. Browsers cannot decode WavPack, so this array holds the DECODED
             samples; the size that WavPack itself needs is recorded in the manifest.
  bits3      the samples with the 3 lowest bits dropped (step 8 counts, about 24 uV). A lossy
             format a browser can decode with nothing but zlib. The array holds the coarse
             values; sample = value * 8 + 3.5.

Every array uses a delta filter along time and zlib, so a chunk is decoded in the browser with
DecompressionStream("deflate") and a running sum. Chunks are (384 channels, 6000 samples).

Also written:
  manifest.json   geometry, scaling, dead channels, variants and their sizes and errors,
                  the high-pass filter coefficients the viewer uses, checksums
  clock_p<probe>.f64   pairs (probe sample, time on the viewer's clock) as float64, see clock_pairs
  test/           reference output of the display processing for one tile (process_tile), which
                  the viewer's ?rawtest=1 compares with its own JavaScript version

Time. The viewer's clock is the behavioral rig's (Timeline), on which the sorted spike times
were placed by piecewise-linear interpolation between sync edges. Timeline samples at 2 kHz, so
that map wobbles by up to about 0.3 ms around a straight line. To draw each spike exactly on
its waveform the viewer needs the same map, and clock_pairs recovers it from the spikes
themselves: every sorted spike has a sample number (Kilosort) and a time (ALF).
"""
import argparse
import hashlib
import json
import time
import zlib
from pathlib import Path

import numpy as np

from cbin import Cbin

FS = 30000
TILE = 6000            # samples per chunk (0.2 s)
MARGIN = 300           # samples of neighbors used when filtering a tile (10 ms)
HP_HZ, HP_ORDER = 300.0, 3
PROBES = ["p1_imec0", "p1_imec1", "p1_imec2", "p1_imec3", "p2_imec0", "p2_imec1", "p2_imec2"]
ALF = Path("D:/temp/fd016_share/alf")
UNITS_NPZ = Path("D:/temp/fd016_raster/units.npz")


def raw_path(session, probe, shank):
    run, imec = PROBES[probe].split("_")            # "p1", "imec0"
    return (Path(session) / "raw_ephys_data" / f"{run}_g0_{imec}" /
            f"{run}_g0_t0.{imec}.sh{shank}.ap.cbin")


def sort_dir(session, probe, shank):
    run, imec = PROBES[probe].split("_")
    return Path(session) / "sorting" / f"{run}_g0_{imec}_shank{shank}"


def clock_pairs(session, probe):
    """(sample, time) for every sorted spike of the probe's QC units, sorted by sample.

    ALF holds each QC unit's spikes inside the synced period, in order; they are a contiguous
    run of that unit's Kilosort spikes, found with a rough linear map and checked at every spike.
    """
    z = np.load(UNITS_NPZ, allow_pickle=True)
    uid = z["unit_ids"]
    assert list(z["probes"]) == PROBES
    fit = json.loads(str(z["sync_fits"]))[probe]
    assert fit["probe"] == PROBES[probe]
    st = np.load(ALF / "spikes.times.npy")
    sc = np.load(ALF / "spikes.clusters.npy")
    mine = np.flatnonzero(uid[:, 0] == probe)
    keep = (sc >= mine.min()) & (sc <= mine.max())
    st, sc = st[keep], sc[keep]
    o = np.argsort(sc, kind="stable")
    st, sc = st[o], sc[o]
    b = np.searchsorted(sc, np.arange(mine.max() + 2))
    PS, PT, skipped = [], [], 0
    for sh in range(4):
        so = sort_dir(session, probe, sh) / "kilosort4" / "sorter_output"
        kt = np.load(so / "spike_times.npy")
        kc = np.load(so / "spike_clusters.npy")
        oo = np.argsort(kc, kind="stable")
        kt, kc = kt[oo], kc[oo]
        kb = np.searchsorted(kc, np.arange(kc.max() + 2))
        for u in np.flatnonzero((uid[:, 0] == probe) & (uid[:, 1] == sh)):
            t = st[b[u]:b[u + 1]]
            s = kt[kb[uid[u, 2]]:kb[uid[u, 2] + 1]]
            if t.size == 0 or s.size == 0:
                skipped += 1
                continue
            g = s / FS * fit["slope"] + fit["offset"]         # rough time of each Kilosort spike
            i0 = np.searchsorted(g, t[0])                     # the Kilosort spike nearest to t[0]
            i0 = min(max(i0 - 1, 0), g.size - 1) if i0 == g.size or (i0 > 0 and t[0] - g[i0 - 1] < g[i0] - t[0]) else i0
            ss = s[i0:i0 + t.size]
            # the rough map is within about 0.6 ms of the true one, so 2 ms at every spike
            # rejects a run that is off by one spike or has a spike missing
            if ss.size != t.size or np.abs(g[i0:i0 + t.size] - t).max() > 0.002:
                skipped += 1
                continue
            PS.append(ss)
            PT.append(t)
    ps = np.concatenate(PS).astype(np.float64)
    pt = np.concatenate(PT)
    o = np.argsort(ps, kind="stable")
    ps, pt = ps[o], pt[o]
    assert np.all(np.diff(pt) >= 0), "time is not monotone in sample number"
    return ps, pt, skipped


def thin_pairs(ps, pt, min_gap=30):
    """Keep pairs at least min_gap samples apart.

    Between sync edges the map is a straight line whose slope differs from 1 by up to a few
    percent, so interpolating across a kept pair's gap is off by a small fraction of a sample
    (measured on this stretch: every spike maps back to within 0.43 sample of its Kilosort
    sample, 99.9% within 0.04).
    """
    keep = [0]
    for i in range(1, ps.size):
        if ps[i] - ps[keep[-1]] >= min_gap:
            keep.append(i)
    return ps[keep], pt[keep]


def highpass_sos():
    from scipy import signal
    return signal.butter(HP_ORDER, HP_HZ, "hp", fs=FS, output="sos")


def process_tile(raw, k, dead, uv_per_count, scale=1.0, offset=0.0):
    """Display processing of tile k of a (channels x samples) array; returns float32 microvolts.

    The viewer's JavaScript (processTile) does the same thing with its own code:
      1. take the tile with up to MARGIN samples on each side; counts = raw * scale + offset
      2. per channel: subtract the mean of that stretch, then high-pass (Butterworth, order
         HP_ORDER at HP_HZ) forward and backward, both from rest
      3. dead channels are set to 0; per sample, subtract the median over the other channels
      4. cut the margins off, multiply by uv_per_count
    """
    from scipy import signal
    n = raw.shape[1]
    a, b = max(0, k * TILE - MARGIN), min(n, (k + 1) * TILE + MARGIN)
    x = np.asarray(raw[:, a:b], np.float64) * scale + offset
    x -= x.mean(axis=1, keepdims=True)
    sos = highpass_sos()
    y = signal.sosfilt(sos, x, axis=1)
    y = signal.sosfilt(sos, y[:, ::-1], axis=1)[:, ::-1]
    good = np.setdiff1d(np.arange(raw.shape[0]), dead)
    y[dead] = 0.0
    y[good] -= np.median(y[good], axis=0, keepdims=True)
    lo = k * TILE - a
    return (y[:, lo:lo + TILE] * uv_per_count).astype(np.float32)


def write_zarr(path, arr):
    """Zarr v2 array, chunks (channels, TILE), delta along time + zlib. Returns bytes stored."""
    import numcodecs
    import zarr
    z = zarr.create_array(str(path), shape=arr.shape, chunks=(arr.shape[0], TILE), dtype="<i2",
                          zarr_format=2, filters=[numcodecs.Delta(dtype="<i2")],
                          compressors=numcodecs.Zlib(level=4), fill_value=0, overwrite=True,
                          chunk_key_encoding={"name": "v2", "separator": "."})
    z[:] = arr
    back = zarr.open_array(str(path), mode="r")[:]
    assert np.array_equal(back, arr), "Zarr round trip changed the data"
    return sum(f.stat().st_size for f in Path(path).iterdir() if f.is_file())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--session", default="Y:/Subjects/FD_016/2026-11-13/1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--probe", type=int, default=0)
    ap.add_argument("--t0", type=float, default=1100.0, help="start, viewer clock (s)")
    ap.add_argument("--dur", type=float, default=60.0)
    ap.add_argument("--shanks", type=int, nargs="*", default=[0, 1, 2, 3])
    a = ap.parse_args()
    from scipy import signal
    from wavpack_numcodecs.wavpack import WavPack
    import wavpack_numcodecs

    out = Path(a.out)
    (out / "test").mkdir(parents=True, exist_ok=True)
    t_begin = time.time()

    # ---- clock: which samples to take, and the map the viewer will use
    ps, pt, skipped = clock_pairs(a.session, a.probe)
    print(f"clock: {ps.size} spike pairs, {skipped} units skipped")
    s_lo = np.interp(a.t0, pt, ps)
    s_hi = np.interp(a.t0 + a.dur, pt, ps)
    c0 = int(s_lo // FS) - 1                      # whole 1 s chunks of the .cbin, 1 s to spare
    c1 = int(s_hi // FS) + 2
    S0, S1 = c0 * FS, c1 * FS
    assert (S1 - S0) % TILE == 0
    m = (ps >= S0 - FS) & (ps <= S1 + FS)
    aps, apt = thin_pairs(ps[m], pt[m])
    np.stack([aps, apt], axis=1).astype("<f8").tofile(out / f"clock_p{a.probe}.f64")
    lin = np.polyfit(aps, apt, 1)
    wob = (apt - np.polyval(lin, aps)) * 1e3
    print(f"samples {S0} to {S1} ({(S1 - S0) / FS:.0f} s); viewer time {np.interp(S0, ps, pt):.3f} to "
          f"{np.interp(S1, ps, pt):.3f} s; {aps.size} clock anchors; deviation from a straight "
          f"line {wob.min():.3f} to {wob.max():.3f} ms")

    sos = highpass_sos()
    variants = {
        "lossless": dict(label="Lossless", scale=1.0, offset=0.0),
        "wv2.25": dict(label="Lossy WavPack, 2.25 bits per sample", scale=1.0, offset=0.0),
        "bits3": dict(label="Lossy, 3 lowest bits dropped", scale=8.0, offset=3.5),
    }
    for v in variants.values():
        v.update(stored_bytes=0, codec_bytes=0, err_sum2=0.0, err_hp_sum2=0.0, err_n=0, err_max=0.0)
    wv_lossy = WavPack(level=3, bps=2.25)
    wv_lossless = WavPack(level=3, bps=None)
    wv_lossless_bytes = 0
    raw_bytes = 0
    shanks = []
    test = None
    for sh in a.shanks:
        cb = Cbin(raw_path(a.session, a.probe, sh))
        assert cb.n_channels == 385 and cb.fs == FS and np.array_equal(cb.bounds[c0:c1 + 1], np.arange(c0, c1 + 1) * FS)
        so = sort_dir(a.session, a.probe, sh)
        prov = json.load(open(so / "provenance.json"))
        pos = np.load(so / "kilosort4" / "sorter_output" / "channel_positions.npy")
        cmap = np.load(so / "kilosort4" / "sorter_output" / "channel_map.npy")
        assert np.array_equal(cmap, np.arange(384)), "sorter dropped or reordered channels"
        dead = sorted(int(i) for i in prov["destripe"]["dead_channel_indices"])
        uv = float(prov["recording"]["uv_per_bit"])
        arrs = {k: np.zeros((384, S1 - S0), np.int16) for k in variants}
        for c in range(c0, c1):
            blk = cb.chunk(c)[:, :384]                       # samples x channels
            sl = slice((c - c0) * FS, (c - c0 + 1) * FS)
            arrs["lossless"][:, sl] = blk.T
            enc = wv_lossy.encode(np.ascontiguousarray(blk))
            dec = wv_lossy.decode(enc)
            dec = (np.frombuffer(dec, "<i2") if isinstance(dec, (bytes, bytearray, memoryview))
                   else np.asarray(dec)).reshape(blk.shape)
            arrs["wv2.25"][:, sl] = dec.T
            variants["wv2.25"]["codec_bytes"] += len(enc)
            wv_lossless_bytes += len(wv_lossless.encode(np.ascontiguousarray(blk)))
            arrs["bits3"][:, sl] = (blk >> 3).T
            raw_bytes += blk.nbytes
            if (c - c0) % 10 == 0:
                print(f"  shank {sh}: read {c - c0 + 1} of {c1 - c0} s ({time.time() - t_begin:.0f} s)")
        found_dead = np.flatnonzero(arrs["lossless"].std(axis=1) == 0).tolist()
        print(f"  shank {sh}: dead per sorter {dead}, zero-variance here {found_dead}")
        # errors of the lossy variants: all samples, and after the display high-pass (10 s).
        # Dead sites are included (6 of 1,536 here).
        ref = arrs["lossless"].astype(np.float32)
        n10 = min(ref.shape[1], 10 * FS)
        ref_hp = signal.sosfiltfilt(sos, ref[:, :n10].astype(np.float64), axis=1)
        for k, v in variants.items():
            if k == "lossless":
                continue
            rec = arrs[k].astype(np.float32) * v["scale"] + v["offset"]
            e = rec - ref
            v["err_sum2"] += float((e.astype(np.float64) ** 2).sum()) * uv * uv
            v["err_max"] = max(v["err_max"], float(np.abs(e).max()) * uv)
            eh = signal.sosfiltfilt(sos, rec[:, :n10].astype(np.float64), axis=1) - ref_hp
            v["err_hp_sum2"] += float((eh ** 2).sum()) * uv * uv
            v["err_n"] += e.size
            v.setdefault("err_hp_n", 0)
            v["err_hp_n"] += eh.size
        noise = None
        y = process_tile(arrs["lossless"], 3, dead, uv)
        noise = float(np.median(np.median(np.abs(y), axis=1)) / 0.6745)
        sums = {}
        for k, v in variants.items():
            v["stored_bytes"] += write_zarr(out / f"p{a.probe}" / f"sh{sh}" / k, arrs[k])
            sums[k] = int(arrs[k][:, :TILE].astype(np.int64).sum())      # checksum of chunk 0
        shanks.append(dict(shank=sh, x_um=pos[:, 0].tolist(), y_um=pos[:, 1].tolist(), dead=dead,
                           uv_per_count=uv, noise_uv=noise, chunk0_sum=sums))
        if test is None:                                   # reference tiles for ?rawtest=1
            kt = 7
            test = dict(shank=sh, tile=kt, files={})
            for k, v in variants.items():
                ref_t = process_tile(arrs[k], kt, dead, uv, v["scale"], v["offset"])
                ref_t.tofile(out / "test" / f"ref_sh{sh}_tile{kt}_{k}.f32")
                test["files"][k] = f"test/ref_sh{sh}_tile{kt}_{k}.f32"
        del arrs
        print(f"  shank {sh} written ({time.time() - t_begin:.0f} s)")

    for k, v in variants.items():
        if k != "lossless":
            v["err_sd_uv"] = (v.pop("err_sum2") / v["err_n"]) ** 0.5
            v["err_hp_sd_uv"] = (v.pop("err_hp_sum2") / v.pop("err_hp_n")) ** 0.5
            v["err_max_uv"] = v.pop("err_max")
        else:
            for kk in ("err_sum2", "err_hp_sum2", "err_max"):
                v.pop(kk)
        v.pop("err_n")
        v["raw_bytes"] = raw_bytes
        if not v["codec_bytes"]:
            v["codec_bytes"] = v["stored_bytes"]            # the array is the real format
        v["ratio"] = raw_bytes / v["codec_bytes"]
    man = dict(
        description="Prototype raw data for the raster viewer's underlay; see make_raw_proto.py",
        probe=a.probe, probe_name=PROBES[a.probe], fs=FS, s0=S0, n_samples=S1 - S0, tile=TILE,
        margin=MARGIN, n_channels=384, clock=f"clock_p{a.probe}.f64",
        highpass=dict(hz=HP_HZ, order=HP_ORDER, sos=sos.tolist()),
        path="p{probe}/sh{shank}/{variant}/0.{tile}",
        variants=[dict(key=k, **v) for k, v in variants.items()],
        wavpack_lossless_ratio=raw_bytes / wv_lossless_bytes,
        wavpack_version=str(wavpack_numcodecs.wavpack_version),
        shanks=shanks, test=test)
    json.dump(man, open(out / "manifest.json", "w"), indent=1)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "label"} for k, v in variants.items()}, indent=1))
    print(f"WavPack lossless ratio {man['wavpack_lossless_ratio']:.2f}; done in {time.time() - t_begin:.0f} s")


if __name__ == "__main__":
    main()
