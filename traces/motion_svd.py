"""Top singular components of the motion energy of a cropped video (after Stringer et al. 2019,
Science 364:eaav7893, as implemented in Facemap).

python motion_svd.py <frames_binned.npz> <out.npz> [n_components] [keyframe_period]

Motion energy: the absolute difference between consecutive frames, per (spatially binned) pixel.
The mean motion-energy image is subtracted, and the singular value decomposition of the
(time x pixel) matrix gives temporal components U[:, k] * S[k] with spatial masks V[k].
Signs are arbitrary in an SVD; each component is flipped so that its mask sums to a positive
number (positive = more motion, where the mask is mostly positive).

keyframe_period (frames, 0 = none): H.264 video stores a full image every so many frames
(30 in the lab's eye videos) and only changes in between, so the frame pairs around each full
image differ a little more than the others even when nothing moves: a comb in the motion energy
at that period (here +31% of the resting level in one pair of 30). With a period given, the
median motion energy of each pixel at each phase of the period, minus that pixel's median over
all phases, is subtracted before anything else. The script prints the resting level by phase
before and after.

Output: temporal (n_frames - 1, K) float32, U * S in units of binned-pixel gray levels;
masks (K, h, w); sing (K,) singular values; var_frac (K,) fraction of the total variance of the
mean-subtracted motion energy; mean_me (h, w); total_me (n_frames - 1,) mean motion energy over
the crop. Sample j belongs to the pair of frames j and j + 1.
"""
import sys

import numpy as np
from scipy.sparse.linalg import svds

src, out = sys.argv[1:3]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 10
PERIOD = int(sys.argv[4]) if len(sys.argv) > 4 else 0
z = np.load(src)
f = z["binned"]
n, h, w = f.shape
me = np.abs(np.diff(f.reshape(n, h * w), axis=0))             # (n - 1, P) float32
del f


def by_phase(x, per):
    """Median of x (time,) at each phase of the period, relative to the overall median."""
    return np.median(x[: (x.size // per) * per].reshape(-1, per), axis=0) / np.median(x)


if PERIOD:
    print("resting motion energy by phase, before:", np.round(by_phase(me.mean(axis=1), PERIOD), 3))
    m = (me.shape[0] // PERIOD) * PERIOD
    excess = np.empty((PERIOD, h * w), np.float32)
    for ph in range(PERIOD):                                   # median over the frame pairs at this phase
        excess[ph] = np.median(me[ph:m:PERIOD], axis=0)
    excess -= np.median(excess, axis=0)
    for ph in range(PERIOD):
        me[ph::PERIOD] -= excess[ph]
    print("resting motion energy by phase, after: ", np.round(by_phase(me.mean(axis=1), PERIOD), 3))
mean_me = me.mean(axis=0)
total_me = me.mean(axis=1)
me -= mean_me
total_var = float((me.astype(np.float64) ** 2).sum())
U, S, Vt = svds(me, k=K, which="LM", tol=1e-8)
o = np.argsort(S)[::-1]
U, S, Vt = U[:, o], S[o], Vt[o]
sgn = np.sign(Vt.sum(axis=1))
sgn[sgn == 0] = 1
U, Vt = U * sgn, Vt * sgn[:, None]
# checks: orthonormal factors, and each triplet satisfies the SVD equations
assert np.allclose(U.T @ U, np.eye(K), atol=1e-4) and np.allclose(Vt @ Vt.T, np.eye(K), atol=1e-4)
resid = np.abs(me @ Vt.T - U * S).max() / S[0]
print("largest residual of M v = s u, relative to s1:", resid)
assert resid < 1e-4
print("singular values:", np.round(S, 1))
print("fraction of variance:", np.round(S ** 2 / total_var, 4), "sum of first 3:", round(float((S[:3] ** 2).sum() / total_var), 4))
np.savez(out, temporal=(U * S).astype(np.float32), masks=Vt.reshape(K, h, w).astype(np.float32), sing=S,
         var_frac=S ** 2 / total_var, mean_me=mean_me.reshape(h, w), total_me=total_me,
         crop=z["crop"], bin=z["bin"], keyframe_period=PERIOD)
