"""Write behavioral traces for the raster viewer (traces.js draws them under the raster).

A trace is a regularly sampled series on the viewer's clock (the clock of the spike times in
data/meta.json): sample i belongs to time t0 + i * dt seconds. Missing data are NaN.

    from prep_traces import write_traces
    write_traces("../data", [dict(name="wheel", label="Wheel velocity", unit="cm/s",
                                  t0=0.005, dt=0.01, values=v)])

writes data/traces.json and one little-endian float32 file per trace:

    {"format": "rasterViewer-traces-1", "traces": [
       {"name": ..., "label": ..., "unit": ..., "file": "trace_<name>.f32", "n": ..., "t0": ...,
        "dt": ..., "ylim": [lo, hi], "n_valid": ..., "sum": ..., "description": ...}, ...]}

name: letters, digits and _ (used in the URL). ylim: the fixed vertical range of the trace's
lane; values beyond it are clipped in the drawing. If not given: the 0.5th to 99.5th percentile,
widened to include zero; for a signal that sits at one value nearly all the time (a wheel that
rarely turns), the percentiles of the samples away from that value.
n_valid and sum (of the float32 values, in float64) let the viewer's self-test (?tracetest=1)
check that it read the same numbers.

python prep_traces.py <data folder>      checks the files of an existing traces.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

FORMAT = "rasterViewer-traces-1"


def default_ylim(v):
    """A vertical range for a lane: robust to rare extremes, always containing zero."""
    x = v[np.isfinite(v)]
    if x.size == 0:
        return [-1.0, 1.0]
    lo, hi = np.percentile(x, [0.5, 99.5])
    if hi <= lo:                                      # at one value at least 99% of the time:
        away = x[x != np.median(x)]                   # use the samples that are not at it
        if away.size:
            lo, hi = np.percentile(away, [0.5, 99.5])
    lo, hi = min(float(lo), 0.0), max(float(hi), 0.0)
    if hi <= lo:
        return [-1.0, 1.0]
    step = 10.0 ** np.floor(np.log10(hi - lo)) / 10    # round outward to two digits of the range
    return [float(np.floor(lo / step) * step), float(np.ceil(hi / step) * step)]


def write_traces(out_dir, traces, note=""):
    out_dir = Path(out_dir)
    entries, seen = [], set()
    for tr in traces:
        name = tr["name"]
        assert re.fullmatch(r"[A-Za-z0-9_]+", name) and name not in seen, f"bad or repeated name {name!r}"
        seen.add(name)
        v = np.asarray(tr["values"], dtype="<f4").ravel()
        assert not np.any(np.isinf(v)), f"{name}: infinite values (use NaN for missing data)"
        assert tr["dt"] > 0 and np.isfinite(tr["t0"])
        ylim = [float(a) for a in tr["ylim"]] if tr.get("ylim") else default_ylim(v)
        assert ylim[1] > ylim[0]
        f = f"trace_{name}.f32"
        v.tofile(out_dir / f)
        ok = np.isfinite(v)
        entries.append(dict(name=name, label=tr.get("label", name), unit=tr.get("unit", ""), file=f,
                            n=int(v.size), t0=float(tr["t0"]), dt=float(tr["dt"]), ylim=ylim,
                            n_valid=int(ok.sum()), sum=float(v[ok].astype(np.float64).sum()),
                            description=tr.get("description", "")))
    man = dict(format=FORMAT, note=note, traces=entries)
    json.dump(man, open(out_dir / "traces.json", "w"), indent=1)
    check(out_dir)
    return man


def check(data_dir):
    """Read traces.json and its files back; assert they agree."""
    data_dir = Path(data_dir)
    man = json.load(open(data_dir / "traces.json"))
    assert man["format"] == FORMAT
    for e in man["traces"]:
        v = np.fromfile(data_dir / e["file"], "<f4")
        ok = np.isfinite(v)
        assert v.size == e["n"] and int(ok.sum()) == e["n_valid"], e["name"]
        assert np.isclose(v[ok].astype(np.float64).sum(), e["sum"], rtol=1e-12, atol=1e-9), e["name"]
        t_end = e["t0"] + (e["n"] - 1) * e["dt"]
        print(f"{e['name']:12s} {e['label']:22s} {e['n']:8d} samples ({e['n_valid']} with data), "
              f"{e['t0']:.4f} to {t_end:.4f} s every {e['dt'] * 1e3:.4f} ms, lane {e['ylim'][0]:g} to {e['ylim'][1]:g} {e['unit']}")
    return man


if __name__ == "__main__":
    check(sys.argv[1])
