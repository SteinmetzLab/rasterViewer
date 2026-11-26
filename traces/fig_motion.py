"""Figure: what the video motion components are (crop, masks, traces, variance).

python fig_motion.py <frames_binned.npz> <motion_svd.npz> <motion_svd_traces.npz> <out.png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

frames_file, svd_file, tr_file, out = sys.argv[1:5]
fz = np.load(frames_file)
mean_frame = fz["binned"][::50].mean(axis=0)
z = np.load(svd_file)
tr = np.load(tr_file)
B = int(z["bin"])
h, w = mean_frame.shape
ext = [0, w * B, h * B, 0]
t = float(tr["t0"]) + np.arange(tr["z"].shape[0]) * float(tr["dt"])

fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(3, 5, height_ratios=[1.15, 1, 1], hspace=0.5, wspace=0.45)
cols = ["#0072B2", "#D55E00", "#009E73"]
OFF, CLIP = 16, 7          # vertical offset between components; values drawn clipped at +-CLIP


def stacked(ax, sl, lw):
    for k in range(3):
        ax.plot(t[sl], np.clip(tr["z"][sl, k], -CLIP, CLIP) - OFF * k, color=cols[k], lw=lw)
        ax.text(0.005, (CLIP + 3 - OFF * k + 2 * OFF + CLIP + 1) / (2 * OFF + 2 * CLIP + 5), f"Component {k + 1}",
                color=cols[k], transform=ax.transAxes, va="top", fontsize=10)
    ax.set_xlim(t[sl][0], t[sl][-1])
    ax.set_ylim(-2 * OFF - CLIP - 1, CLIP + 4)
    ax.set_yticks([-2 * OFF - 5, -2 * OFF, -2 * OFF + 5, -OFF - 5, -OFF, -OFF + 5, -5, 0, 5])
    ax.set_yticklabels(["−5", "0", "5"] * 3)
    ax.set_title("")
    ax.set_xlabel("Time in the session (s)")
    ax.set_ylabel(f"Motion component (sd from rest,\nclipped at ±{CLIP})")

ims = [("Mean frame in the crop", mean_frame, "gray", None),
       ("Mean motion energy", z["mean_me"], "magma", None)] + \
      [(f"Mask of component {k + 1}", z["masks"][k], "RdBu_r", "sym") for k in range(3)]
for i, (title, im, cmap, mode) in enumerate(ims):
    ax = fig.add_subplot(gs[0, i])
    kw = dict(vmin=-np.abs(im).max(), vmax=np.abs(im).max()) if mode else {}
    ax.imshow(im, cmap=cmap, extent=ext, **kw)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x in the crop (pixels)")
    if i == 0:
        ax.set_ylabel("y in the crop (pixels)")
    else:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
ax = fig.add_subplot(gs[1, :4])
stacked(ax, slice(None), 0.5)
ax.set_title("Time course of the three components, whole video", fontsize=11)
ax.text(0.995, 0.98, f"{tr['z'].shape[0]:,} frame pairs", transform=ax.transAxes, ha="right", va="top", fontsize=10)
ax = fig.add_subplot(gs[1, 4])
vf = z["var_frac"]
ax.plot(np.arange(1, vf.size + 1), 100 * vf, "o-", color="k", ms=4)
ax.set_yscale("log")
ax.set_xlabel("Component")
ax.set_ylabel("Variance of motion energy (%)")
ax.set_title("Variance per component", fontsize=11)
ax.text(0.97, 0.97, f"First 3: {100 * vf[:3].sum():.1f}%", transform=ax.transAxes, ha="right", va="top", fontsize=10)
ax = fig.add_subplot(gs[2, :4])
i0 = int(np.argmax(np.convolve((np.abs(tr["z"]) > 1).sum(axis=1), np.ones(700), mode="valid")))   # busiest 10 s
stacked(ax, slice(i0, i0 + 700), 0.9)
ax.set_title("The 10 s with the most motion", fontsize=11)
ax = fig.add_subplot(gs[2, 4])
c = np.corrcoef(np.column_stack([tr["z"], np.log(z["total_me"] + 1e-3)]).T)
ax.imshow(c, cmap="RdBu_r", vmin=-1, vmax=1)
lab = ["1", "2", "3", "Total"]
ax.set_xticks(range(4)); ax.set_xticklabels(lab); ax.set_yticks(range(4)); ax.set_yticklabels(lab)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{c[i, j]:.2f}", ha="center", va="center", fontsize=9, color="k")
ax.set_xlabel("Component (Total: log mean motion energy)")
ax.set_ylabel("Component")
ax.set_title("Correlation between the traces", fontsize=11)
fig.savefig(out, dpi=130, bbox_inches="tight")
print("saved", out)
