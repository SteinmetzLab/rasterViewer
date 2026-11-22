"""Copy the Rastermap ordering into the viewer's data folder (run once, instant).

Input:  rastermap/rastermap_main.json (written by rastermap/run_rastermap.py; see the
        README there). isort[k] is the viewer display row shown at position k.
Output: data/rastermap.json with isort, the fit parameters, and the source file name.

Checks that isort is a permutation of the viewer's rows and that the cluster ids match.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "rastermap" / "rastermap_main.json"
OUT = HERE / "data" / "rastermap.json"


def main():
    rm = json.load(open(SRC))
    meta = json.load(open(HERE / "data" / "meta.json"))
    n = meta["n_rows"]
    isort = [int(x) for x in rm["isort"]]
    assert sorted(isort) == list(range(n)), "isort is not a permutation of the viewer's rows"
    assert rm["cluster_id"] == meta["rows"]["cluster"], "cluster ids differ from meta.json"
    json.dump(dict(source=f"rastermap/{SRC.name}", params=rm["params"], isort=isort),
              open(OUT, "w"))
    print(f"wrote {OUT} ({n} rows)")


if __name__ == "__main__":
    main()
