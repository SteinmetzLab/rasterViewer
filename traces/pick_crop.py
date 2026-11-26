"""Pop up frames of a video and let the user draw the crop rectangle.

python pick_crop.py <sample_frames.npy> <crop.json>

sample_frames.npy: (n, height, width) uint8, a few frames spread over the video.
Writes crop.json {"x0", "y0", "x1", "y1"} in pixels (x1 and y1 exclusive) when OK is clicked.
Exit code 0 after OK, 1 if the window was closed without OK.
"""
import json
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, RectangleSelector

plt.rcParams["font.family"] = "Arial"
frames = np.load(sys.argv[1])
out = sys.argv[2]
n, H, W = frames.shape
state = dict(i=n // 2, box=None, ok=False)

fig, ax = plt.subplots(figsize=(11, 9))
fig.subplots_adjust(bottom=0.12, top=0.93)
im = ax.imshow(frames[state["i"]], cmap="gray", vmin=0, vmax=255)
ax.set_xlabel("x (pixels)")
ax.set_ylabel("y (pixels)")
info = fig.text(0.5, 0.955, "", ha="center", fontsize=12)


def refresh():
    b = state["box"]
    s = "Drag a rectangle on the image (you can redraw it or drag its handles), then click OK."
    if b:
        s = (f"Crop: x {b[0]} to {b[2]}, y {b[1]} to {b[3]}  ({b[2] - b[0]} x {b[3] - b[1]} pixels). "
             "Adjust if needed, then click OK.")
    info.set_text(s + f"   Showing sample frame {state['i'] + 1} of {n}.")
    fig.canvas.draw_idle()


def on_select(e0, e1):
    x0, x1 = sorted((e0.xdata, e1.xdata))
    y0, y1 = sorted((e0.ydata, e1.ydata))
    x0, y0 = max(0, int(round(x0))), max(0, int(round(y0)))
    x1, y1 = min(W, int(round(x1))), min(H, int(round(y1)))
    state["box"] = (x0, y0, x1, y1) if (x1 - x0 >= 8 and y1 - y0 >= 8) else None
    refresh()


sel = RectangleSelector(ax, on_select, useblit=False, interactive=True, button=[1],
                        props=dict(edgecolor="#00E0FF", facecolor="none", linewidth=2))


def on_next(_):
    state["i"] = (state["i"] + 1) % n
    im.set_data(frames[state["i"]])
    refresh()


def on_ok(_):
    if not state["box"]:
        info.set_text("No rectangle yet: drag one on the image first.")
        fig.canvas.draw_idle()
        return
    x0, y0, x1, y1 = state["box"]
    json.dump(dict(x0=x0, y0=y0, x1=x1, y1=y1, frame_width=W, frame_height=H), open(out, "w"))
    state["ok"] = True
    plt.close(fig)


b_next = Button(fig.add_axes([0.30, 0.02, 0.18, 0.05]), "Show another frame")
b_ok = Button(fig.add_axes([0.52, 0.02, 0.18, 0.05]), "OK")
b_next.on_clicked(on_next)
b_ok.on_clicked(on_ok)
refresh()
try:
    fig.canvas.manager.set_window_title("FD_016 eye video: choose the crop for motion energy")
    win = fig.canvas.manager.window
    if matplotlib.get_backend().lower().startswith("qt"):
        win.raise_(); win.activateWindow()
    else:
        win.attributes("-topmost", True)
except Exception:
    pass
plt.show()
sys.exit(0 if state["ok"] else 1)
