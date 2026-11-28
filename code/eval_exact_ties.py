"""Exact expected values of the top-K evidence measures under random tie-breaking.

Many drugs share a score, so which of them fall inside the top K is random. eval_metrics.py and eval_tables.py estimate
the average by sampling 500 tie orders; close to a rounding boundary the sampled value can round differently. This
script computes the average over all tie orders directly. The ranking-quality, ablation and grade-sensitivity values
reported in the paper come from this script. The sampled result files are kept unchanged.

For a group of n tied drugs with s places left inside the top K, each drug is included with probability s/n, and each
of the s positions holds, on average, the mean gain of the group.

SP is left out: 233 of 299 drugs share its top score. The evidence grades are AI-assisted and have not been verified
by a human expert.

    python code/eval_exact_ties.py     # writes results/Evaluation/main/exact_tie_values.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from evalcore import Data, adjacency, pathway_links  # noqa: E402
from eval_run import all_methods, ranked  # noqa: E402
from eval_tables import ALIAS, HARM, parent  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "results/Evaluation/evidence"
OUT = ROOT / "results/Evaluation/main/exact_tie_values.csv"

METHODS = ["NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed", "Decidable-only", "Consistent-any", "Direction-aware DWPC"]
BORDERLINE = ["CLOFIBRATE", "MITOXANTRONE", "MIBAVADEMAB", "LABETALOL"]
CONFLICTING = ["CELECOXIB", "CIPROFIBRATE", "CILENGITIDE", "PACLITAXEL"]
SCENARIOS = {
    "base": {},
    "borderline four set to 0": {p: 0 for p in BORDERLINE},
    "borderline four set to 1": {p: 1 for p in BORDERLINE},
    "borderline and conflicting eight set to 0": {p: 0 for p in BORDERLINE + CONFLICTING},
}


def exact(table, k, grades, ideal):
    """Expected measures of the top k of `table` (columns score, parent) over all orders of tied scores."""
    t = table.assign(grade=table.parent.map(grades).fillna(0)).sort_values("score", ascending=False)
    top_ideal = np.asarray(ideal[:k], dtype=float)
    idcg = ((2 ** top_ideal - 1) / np.log2(np.arange(len(top_ideal)) + 2)).sum()
    size = min(k, len(t))
    out = dict.fromkeys(["support", "ndcg", "human", "animal", "in_vitro", "harm"], 0.0)
    pos = 0
    for _, group in t.groupby("score", sort=False):
        slots = min(len(group), max(k - pos, 0))
        if slots == 0:
            break
        w = slots / len(group)
        out["support"] += w * (group.grade >= 1).sum() / size
        out["human"] += w * (group.grade == 3).sum()
        out["animal"] += w * (group.grade == 2).sum()
        out["in_vitro"] += w * (group.grade == 1).sum()
        out["harm"] += w * group.parent.isin(HARM).sum()
        out["ndcg"] += (2 ** group.grade - 1).mean() * (1 / np.log2(np.arange(pos, pos + slots) + 2)).sum() / idcg
        pos += slots
    return out


def main():
    D = Data()
    scores = all_methods(D, adjacency(D), pathway_links(D, D.seeds), D.seeds.i.to_numpy())[0]
    base = pd.read_csv(EV / "evidence_grades.csv").set_index("parent").grade.to_dict()
    pool = pd.read_csv(EV / "judged_pool.csv")
    rows = []
    for label, changes in SCENARIOS.items():
        grades = {**base, **changes}
        ideal = sorted((changes.get(p, g) for p, g in zip(pool.parent, pool.grade)), reverse=True)
        for m in METHODS:
            t = ranked(scores[m])
            t = t.assign(parent=t.Drug_ID.map(D.name).map(parent).replace(ALIAS))
            for k in (10, 20):
                rows.append({"scenario": label, "method": m, "K": k, **exact(t, k, grades, ideal), "ranked": len(t)})
    out = pd.DataFrame(rows)
    out.round(4).to_csv(OUT, index=False)
    print(out[(out.scenario == "base") & (out.K == 20)].round(4).to_string(index=False))
    print("written", OUT)


if __name__ == "__main__":
    main()
