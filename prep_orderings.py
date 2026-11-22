"""Copy the further neuron orderings into the viewer's data folder (run once, instant).

Input:  rastermap/orderings.json (written by rastermap/make_orderings.py; see the README there).
Output: data/orderings.json: a list of {key, label, axis, description, params, isort}. The
        viewer adds one entry to its Neurons menu per ordering (layouts of the same kind as
        "Rastermap order": the neuron at position k is display row isort[k]).

Checks that every isort is a permutation of the viewer's rows and that the cluster ids match.
The plain Rastermap order stays in data/rastermap.json (prep_rastermap.py).
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "rastermap" / "orderings.json"
OUT = HERE / "data" / "orderings.json"

# key in rastermap/orderings.json -> (menu label, axis label after "Neuron, ")
MENU = [
    ("rm_refine", "Rastermap, refined", "Rastermap order, locally refined"),
    ("rm_tree", "Tree (Rastermap layout)", "clustering-tree order (branches follow Rastermap)"),
    ("hp_rastermap", "Rastermap, no slow drift",
     "Rastermap order (sorted on periods under {p:g} s only)"),
    ("hp_refine", "Rastermap, refined, no slow drift",
     "Rastermap order, locally refined (sorted on periods under {p:g} s only)"),
    ("hp_tree", "Tree, no slow drift",
     "clustering-tree order (sorted on periods under {p:g} s only)"),
]


def main():
    src = json.load(open(SRC))
    meta = json.load(open(HERE / "data" / "meta.json"))
    n = meta["n_rows"]
    assert src["cluster_id"] == meta["rows"]["cluster"], "cluster ids differ from meta.json"
    out = []
    for key, label, axis in MENU:
        o = src["orderings"][key]
        isort = [int(x) for x in o["isort"]]
        assert sorted(isort) == list(range(n)), f"{key}: isort is not a permutation of the rows"
        p = o["params"].get("highpass_period_s") or 0
        out.append(dict(key=key, label=label, axis=axis.format(p=p), description=o["description"],
                        params=o["params"], isort=isort))
    json.dump(dict(source=f"rastermap/{SRC.name}", orderings=out), open(OUT, "w"))
    print(f"wrote {OUT} ({len(out)} orderings of {n} rows)")


if __name__ == "__main__":
    main()
