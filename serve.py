"""Serve the viewer locally: python serve.py  ->  http://localhost:8765/

Serves this folder: index.html, viewer.js, data/ (written by prep_data.py) and refs/
(written by make_refs.py, for ?test=1). Works on Windows and macOS.
"""
import argparse
import functools
import http.server
from pathlib import Path

HERE = Path(__file__).resolve().parent
MOUNTS = {}   # everything lives under this folder


class Handler(http.server.SimpleHTTPRequestHandler):
    def translate_path(self, path):
        p = path.split("?", 1)[0].split("#", 1)[0]
        for prefix, root in MOUNTS.items():
            if p.startswith(prefix):
                rel = p[len(prefix):]
                target = (root / rel).resolve()
                if root.resolve() in target.parents or target == root.resolve():
                    return str(target)
                return str(root / "__forbidden__")
        return super().translate_path(path)

    def end_headers(self):
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, fmt, *args):
        pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args()
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", a.port),
                                          functools.partial(Handler, directory=str(HERE)))
    print(f"Serving http://localhost:{a.port}/  (Ctrl+C to stop)")
    srv.serve_forever()
