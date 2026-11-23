"""Read pieces of an mtscomp .cbin file (the lab's compressed SpikeGLX format) without mtscomp.

A .cbin is a sequence of independently compressed chunks (here 1 s each); the .ch file is JSON
with the chunk boundaries (in samples) and byte offsets. Each chunk is zlib-compressed int16,
C order (samples x channels), stored as differences along time when do_time_diff is set.
"""
import json
import zlib
from pathlib import Path

import numpy as np


class Cbin:
    def __init__(self, cbin_path):
        self.path = Path(cbin_path)
        ch = json.load(open(self.path.with_suffix(".ch")))
        assert ch["algorithm"] == "zlib" and ch["dtype"] == "int16" and ch["chunk_order"] == "C"
        assert not ch["do_spatial_diff"], "spatial differences are not handled"
        self.time_diff = bool(ch["do_time_diff"])
        self.n_channels = int(ch["n_channels"])
        self.n_samples = int(ch["shape"][0])
        self.fs = float(ch["sample_rate"])
        self.bounds = np.asarray(ch["chunk_bounds"], np.int64)
        self.offsets = np.asarray(ch["chunk_offsets"], np.int64)
        assert self.bounds.size == self.offsets.size and self.bounds[-1] == self.n_samples

    @property
    def n_chunks(self):
        return self.bounds.size - 1

    def chunk(self, k):
        """Samples bounds[k] .. bounds[k+1] as int16 (n, n_channels)."""
        with open(self.path, "rb") as f:
            f.seek(int(self.offsets[k]))
            buf = f.read(int(self.offsets[k + 1] - self.offsets[k]))
        a = np.frombuffer(zlib.decompress(buf), dtype="<i2")
        n = int(self.bounds[k + 1] - self.bounds[k])
        assert a.size == n * self.n_channels, "chunk has the wrong size"
        a = a.reshape(n, self.n_channels)
        if self.time_diff:
            a = np.cumsum(a, axis=0, dtype=np.int16)      # int16 wrap-around undoes np.diff
        return a
