"""Do neuron orders from different k-means seeds agree better when the cluster order is searched
to (near) its best value? Seeds 0-7, first half of the session, 100 ms bins.
greedy = Rastermap default; best = annealing (16 runs) + exhaustive polish (fork's fastsort).
Writes seed_consistency.json (pairwise |Spearman| of neuron positions, and the share of each
neuron's 20 nearest neighbors kept)."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from rastermap import Rastermap, fastsort
from rastermap.cluster import scaled_kmeans
from validate_rastermap import neighbor_overlap, rank_of

HERE = Path(__file__).resolve().parent
z = np.load("D:/temp/fd016_rastermap_cc_seed0.npz")     # PCs of the first half (saved by a tuning run on the lab desktop)
Usv, Vsv = z["Usv"], z["Vsv"]
good = ~np.isnan(Usv[:, 0])
orders = {"greedy": [], "best": []}
obj = {"greedy": [], "best": []}
for seed in range(8):
    m = Rastermap(random_state=seed, verbose=False).fit(Usv=Usv, Vsv=Vsv)
    orders["greedy"].append(np.asarray(m.isort))
    obj["greedy"].append(float((m.cc.astype(np.float64) * m.BBt).sum()))
    U0 = scaled_kmeans(Usv[good], n_clusters=100, random_state=seed)[0]
    cc = U0 @ U0.T
    _, inds, BBt, info = fastsort.anneal_sort(cc, n_restarts=16, random_state=seed, do_polish=True)
    m2 = Rastermap(random_state=seed, verbose=False).fit(Usv=Usv, Vsv=Vsv, U_nodes=U0[inds])
    F2 = float((m2.cc.astype(np.float64) * m2.BBt).sum())
    assert F2 >= info["objective"] - 1e-2, (F2, info["objective"])   # greedy keeps a polished order
    orders["best"].append(np.asarray(m2.isort))
    obj["best"].append(F2)
    print(seed, round(obj["greedy"][-1], 2), round(F2, 2), flush=True)
out = {"cluster_objective": obj}
for k, O in orders.items():
    R = [rank_of(o) for o in O]
    rho, ov = [], []
    for i in range(len(R)):
        for j in range(i + 1, len(R)):
            r = spearmanr(R[i], R[j]).statistic
            rj = R[j] if r > 0 else len(R[j]) - 1 - R[j]
            rho.append(abs(float(r)))
            ov.append(neighbor_overlap(R[i], rj, 20))
    out[k] = dict(spearman_mean=round(float(np.mean(rho)), 3), spearman_min=round(float(np.min(rho)), 3),
                  spearman_max=round(float(np.max(rho)), 3), neighbor20_mean=round(float(np.mean(ov)), 3))
    print(k, out[k], flush=True)
json.dump(out, open(HERE / "seed_consistency.json", "w"), indent=1)
