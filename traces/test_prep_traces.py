"""Tests of prep_traces.py:  python test_prep_traces.py"""
import tempfile

import numpy as np

from prep_traces import check, default_ylim, write_traces

rng = np.random.default_rng(0)
x = rng.standard_normal(10000)
lo, hi = default_ylim(x)
assert -3.2 < lo < -2.3 and 2.3 < hi < 3.2, (lo, hi)
x1 = x.copy(); x1[5] = 1e6                               # one wild value does not set the range
assert default_ylim(x1) == [lo, hi] or abs(default_ylim(x1)[1] - hi) < 0.2, default_ylim(x1)
w = np.zeros(100000); w[500:520] = -8; w[9000:9030] = 6    # a signal at one value nearly always
assert default_ylim(w) == [-8.0, 6.0], default_ylim(w)
assert default_ylim(np.full(50, np.nan)) == [-1.0, 1.0]
assert default_ylim(np.zeros(50)) == [-1.0, 1.0]
a, b = default_ylim(np.abs(x))                              # positive only: the lane starts at zero
assert a == 0.0 and b > 2
a, b = default_ylim(-np.abs(x))
assert b == 0.0 and a < -2
with tempfile.TemporaryDirectory() as d:
    v = x.astype(np.float32); v[100:200] = np.nan
    man = write_traces(d, [dict(name="a1", label="A", unit="z", t0=1.5, dt=0.01, values=v),
                           dict(name="b", t0=0, dt=1, values=[1, 2, 3], ylim=[0, 5])])
    assert man["traces"][0]["n_valid"] == 9900 and man["traces"][1]["ylim"] == [0.0, 5.0]
    back = np.fromfile(d + "/trace_a1.f32", "<f4")
    assert np.array_equal(back, v, equal_nan=True)
    check(d)
    for bad in (dict(name="a b", t0=0, dt=1, values=[1]), dict(name="c", t0=0, dt=0, values=[1]),
                dict(name="d", t0=0, dt=1, values=[np.inf])):
        try:
            write_traces(d, [bad]); raise SystemExit("accepted a bad trace: %r" % bad)
        except AssertionError:
            pass
print("prep_traces tests passed")
