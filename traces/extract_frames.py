"""Read every frame of a video once; save the cropped, spatially binned frames and coarse
whole-frame brightness (for checking the timing against Timeline signals).

python extract_frames.py <video> <crop.json> <out.npz> [bin]

Output: binned (n_frames, h, w) float32, mean of bin x bin pixel blocks inside the crop (the crop
is trimmed at its right and bottom edges to a multiple of bin); grid (n_frames, 8, 8) float32,
mean brightness of the whole frame in an 8 x 8 grid of blocks; crop (x0, y0, x1, y1) as used.
"""
import json
import sys

import cv2
import numpy as np

video, crop_file, out = sys.argv[1:4]
B = int(sys.argv[4]) if len(sys.argv) > 4 else 4
c = json.load(open(crop_file))
x0, y0 = c["x0"], c["y0"]
w, h = (c["x1"] - x0) // B, (c["y1"] - y0) // B
x1, y1 = x0 + w * B, y0 + h * B
v = cv2.VideoCapture(video)
n = int(v.get(cv2.CAP_PROP_FRAME_COUNT))
H, W = int(v.get(cv2.CAP_PROP_FRAME_HEIGHT)), int(v.get(cv2.CAP_PROP_FRAME_WIDTH))
binned = np.zeros((n, h, w), np.float32)
grid = np.zeros((n, 8, 8), np.float32)
i = 0
while True:
    ok, im = v.read()
    if not ok:
        break
    assert i < n, "more frames than the header says"
    g = im[..., 0]                                   # gray video: the three planes are equal
    binned[i] = g[y0:y1, x0:x1].reshape(h, B, w, B).mean(axis=(1, 3))
    grid[i] = g.reshape(8, H // 8, 8, W // 8).mean(axis=(1, 3))
    i += 1
    if i % 2000 == 0:
        print(i, flush=True)
print("frames read", i, "header says", n)
np.savez(out, binned=binned[:i], grid=grid[:i], crop=np.array([x0, y0, x1, y1]), bin=B)
