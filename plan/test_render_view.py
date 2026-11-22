"""Checks of render_view.Raster.view against brute force (run: python test_render_view.py)."""
import numpy as np
from render_view import Raster

R = Raster()
rng = np.random.default_rng(0)
fails = 0
for trial in range(6):
    t0 = rng.uniform(5, 2300); t1 = t0 + rng.choice([0.05, 3.0, 70.0])
    n0 = rng.uniform(0, 4000); n1 = min(R.n, n0 + rng.choice([7.3, 40.0, 300.5]))
    W, H = int(rng.integers(20, 90)), int(rng.integers(5, 40))
    img = R.view(t0, t1, n0, n1, W, H)
    # brute force: loop over pixels, each spike counted with its row's overlap fraction
    m = (R.t >= t0) & (R.t < t1)
    t, r = R.t[m], R.r[m]
    ref = np.zeros((H, W)); dt = (t1 - t0) / W; ye = np.linspace(n0, n1, H + 1)
    for k in range(H):
        rows = np.arange(int(np.floor(ye[k])), int(np.ceil(ye[k + 1])))
        ov = np.clip(np.minimum(ye[k + 1], rows + 1) - np.maximum(ye[k], rows), 0, None)
        for j in range(W):
            mm = (t >= t0 + j * dt) & (t < t0 + (j + 1) * dt)
            c = np.array([np.sum(r[mm] == q) for q in rows])
            ref[k, j] = (ov * c).sum() / (ov.sum() * dt)
    err = np.max(np.abs(img - ref) / np.maximum(ref, 1e-9))
    # conservation: sum(rate * coverage * dt) over the image = spikes in rows/time, weighted
    ok = np.allclose(img, ref, rtol=1e-6, atol=1e-9)
    fails += not ok
    print(f"t {t0:.1f}+{t1-t0:g}s rows {n0:.1f}-{n1:.1f} {W}x{H}: max rel err {err:.2e} {'OK' if ok else 'FAIL'}")
# whole-image conservation on the full view
img = R.view(R.t[0], R.t[-1] + 1e-9, 0, R.n, 300, 100)
dt = (R.t[-1] + 1e-9 - R.t[0]) / 300
total = img.sum() * dt * (R.n / 100)
print(f"conservation: {total:.0f} vs {len(R.t)} spikes, rel diff {abs(total-len(R.t))/len(R.t):.1e}")
print("ALL OK" if fails == 0 else f"{fails} FAILED")
