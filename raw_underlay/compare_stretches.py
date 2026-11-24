"""The motion shift, tested on the same units at two times.

    python compare_stretches.py A.json B.json

A.json and B.json are outputs of validate_raw.py for the same shank in two stretches, one in
which the shank had moved and one in which it had not. For every unit with a spike-triggered
average of 30 uV or more in both, the waveform's depth in the probe's own coordinates (y_raw)
is compared between the stretches. The unit's depth in the sort (motion-corrected) is the same
in both, so if the viewer's rule is right (a site at y is drawn at y - d),

    y_raw(A) - y_raw(B) = d(A) - d(B)

for every unit, whatever bias the depth estimate has for that unit. Reports the slope of the
left side on the right side (1 if the rule is right, 0 if the motion estimate were unrelated
to the data, -1 for the opposite sign) and what remains under each of the three.
"""
import json
import sys

import numpy as np


def main():
    A, B = (json.load(open(f)) for f in sys.argv[1:3])
    assert (A["summary"]["probe"], A["summary"]["shank"]) == (B["summary"]["probe"], B["summary"]["shank"])
    b = {u["row"]: u for u in B["units"]}
    pairs = [(u, b[u["row"]]) for u in A["units"] if u["row"] in b and u["amp_uv"] >= 30 and b[u["row"]]["amp_uv"] >= 30]
    dy = np.array([u["y_raw_um"] - v["y_raw_um"] for u, v in pairs])
    dd = np.array([u["motion_um"] - v["motion_um"] for u, v in pairs])
    med = lambda e: float(np.median(np.abs(e)))
    out = dict(probe=A["summary"]["probe"], shank=A["summary"]["shank"],
               stretches_s=[A["summary"]["stretch_s"], B["summary"]["stretch_s"]], n_units=len(pairs),
               motion_difference_um=dict(median_abs=med(dd), max_abs=float(np.abs(dd).max())),
               slope=float(np.polyfit(dd, dy, 1)[0]), slope_through_origin=float((dd * dy).sum() / (dd * dd).sum()),
               correlation=float(np.corrcoef(dd, dy)[0, 1]),
               median_abs_remainder_um={"viewer's rule": med(dy - dd), "no shift": med(dy), "opposite sign": med(dy + dd)})
    for thr in (5, 10):
        m = np.abs(dd) >= thr
        if m.sum() >= 3:
            out[f"units_moved_{thr}um"] = dict(
                n=int(m.sum()), slope_through_origin=float((dd[m] * dy[m]).sum() / (dd[m] * dd[m]).sum()),
                median_abs_remainder_um={"viewer's rule": med((dy - dd)[m]), "no shift": med(dy[m]),
                                         "opposite sign": med((dy + dd)[m])})
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
