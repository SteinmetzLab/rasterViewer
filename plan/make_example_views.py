"""Example stills of what the multiscale raster viewer would show at several zoom levels.

Each image is rendered pixel-exact by render_view.Raster.view: one image pixel = one screen
pixel of the axes, value = mean rate over all neurons and time inside it.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from render_view import Raster, T_START, T_END

plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 11})
OUT = Path(__file__).parent / "figs"
OUT.mkdir(exist_ok=True)
R = Raster()
N = R.n
W, H, DPI = 1400, 800, 100

# probe 3 (index 2) shank 0 rows, for the "few neurons" views
rows_p3 = np.flatnonzero((R.probe_of_row == 2) & (R.shank_of_row == 0))
p0 = rows_p3[0]
T0 = 1102.0     # the movie's least-synchronized window

VIEWS = [
    # name, t0, t1, n0, n1, title, normalize
    ("01_full_session_all_neurons", T_START, T_END, 0, N,
     "Whole session, all neurons", False),
    ("02_full_session_all_neurons_normalized", T_START, T_END, 0, N,
     "Whole session, all neurons, each neuron divided by its mean rate", True),
    ("03_60s_all_neurons", T0, T0 + 60, 0, N, "60 s, all neurons", False),
    ("04_4s_all_neurons", T0, T0 + 4, 0, N, "4 s, all neurons", False),
    ("05_full_session_one_shank", T_START, T_END, p0, p0 + len(rows_p3),
     f"Whole session, one shank ({len(rows_p3)} neurons, Probe 3 shank 0)", False),
    ("06_full_session_one_shank_normalized", T_START, T_END, p0, p0 + len(rows_p3),
     f"Whole session, one shank, each neuron divided by its mean rate", True),
    ("07_4s_one_shank", T0, T0 + 4, p0, p0 + len(rows_p3),
     "4 s, one shank", False),
    ("08_500ms_60_neurons", T0, T0 + 0.5, p0, p0 + 60, "500 ms, 60 neurons", False),
    ("09_50ms_20_neurons", T0, T0 + 0.05, p0, p0 + 20, "50 ms, 20 neurons", False),
    ("10_50ms_20_neurons_ticks", T0, T0 + 0.05, p0, p0 + 20,
     "50 ms, 20 neurons, raster mode: each spike drawn at least 5 pixels wide", "ticks"),
]

probe_mid = [np.mean(np.flatnonzero(R.probe_of_row == p)) for p in range(7)]


def fmt_scale(t0, t1, n0, n1):
    dt = (t1 - t0) / W
    dn = (n1 - n0) / H
    tdt = f"{dt:.2f} s" if dt >= 0.1 else (f"{dt*1e3:.2f} ms" if dt >= 1e-4 else f"{dt*1e6:.0f} µs")
    return dt, dn, f"1 pixel = {dn:.2f} neurons × {tdt}"


stats = {}
for name, t0, t1, n0, n1, title, norm in VIEWS:
    ticks = norm == "ticks"
    norm = False if ticks else norm
    img = R.view(t0, t1, n0, n1, W, H, normalize=norm)
    if ticks:   # display-only dilation so 1-pixel spikes stay visible
        from scipy.ndimage import maximum_filter1d
        img = maximum_filter1d(img, 5, axis=1)
    dt, dn, scale = fmt_scale(t0, t1, n0, n1)
    nz = img[img > 0]
    stats[name] = dict(px_dt_s=dt, px_neurons=dn, frac_zero=float((img == 0).mean()),
                       p1=float(np.percentile(nz, 1)) if nz.size else None,
                       p50=float(np.percentile(nz, 50)) if nz.size else None,
                       p99=float(np.percentile(nz, 99)) if nz.size else None,
                       max=float(img.max()))
    fig = plt.figure(figsize=((W + 260) / DPI, (H + 150) / DPI), dpi=DPI)
    ax = fig.add_axes([150 / (W + 260), 80 / (H + 150), W / (W + 260), H / (H + 150)])
    if norm:
        cnorm = LogNorm(0.1, 10); cmap = plt.get_cmap("RdBu_r").copy()
        clab = "Rate / neuron's session mean"
    else:
        cnorm = LogNorm(0.3, 300); cmap = plt.get_cmap("magma").copy()
        clab = "Mean firing rate (spikes/s)"
    cmap.set_bad("black"); cmap.set_under(cmap(0.0) if not norm else "#08306b")
    shown = np.where(img > 0, img, np.nan) if not norm else np.where(img > 0, img, 1e-6)
    long_t = (t1 - t0) > 10
    if long_t:
        ext = [t0, t1, n1, n0]; xl = "Time (s)"
    else:
        ext = [0, (t1 - t0) * 1e3, n1, n0]; xl = f"Time from {t0:.0f} s (ms)"
    im = ax.imshow(shown, aspect="auto", interpolation="nearest", cmap=cmap, norm=cnorm,
                   extent=ext)
    ax.set_xlabel(xl)
    if (n1 - n0) > 1000:
        ax.set_yticks(probe_mid, [f"Probe {p+1}" for p in range(7)])
        ax.set_ylabel("Neuron (probe, shank, depth order)")
    else:
        ax.set_ylabel("Neuron (display row)")
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_title(title, loc="left", fontsize=13)
    ax.text(0.995, 1.012, scale, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=10, color="0.3")
    cax = fig.add_axes([(W + 165) / (W + 260), 80 / (H + 150), 16 / (W + 260), H / (H + 150)])
    cb = fig.colorbar(im, cax=cax, extend="both")
    cb.set_label(clab)
    fig.savefig(OUT / f"{name}.png", dpi=DPI)
    plt.close(fig)
    print(name, json.dumps({k: (round(v, 4) if isinstance(v, float) else v)
                            for k, v in stats[name].items()}))

json.dump(stats, open(OUT / "view_stats.json", "w"), indent=1)
