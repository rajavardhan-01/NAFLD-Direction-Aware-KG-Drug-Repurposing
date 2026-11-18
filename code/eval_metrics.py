"""Evidence support and judged-pool NDCG of every method (Table A) from the graded literature screen.

grade 3 human clinical, 2 animal in vivo, 1 in vitro, 0 no positive support found (unknown, harm signals are flagged
separately in evidence_grades.csv). A drug inherits the grade of its parent compound (salt forms share it).

    python code/eval_metrics.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from eval_evidence_fetch import parent  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "results/Evaluation/evidence"
MAIN = ROOT / "results/Evaluation/main"
METHODS = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]
HARM = {"DAUNORUBICIN", "DOXIFLURIDINE", "FLOXURIDINE", "CERTOLIZUMAB", "GOLIMUMAB"}


def dcg(rel):
    rel = np.asarray(rel, dtype=float)
    return float(((2 ** rel - 1) / np.log2(np.arange(len(rel)) + 2)).sum())


def main():
    r = pd.read_csv(MAIN / "all_rankings.csv")
    g = pd.read_csv(EV / "evidence_grades.csv").set_index("parent").grade.to_dict()
    r["parent"] = r.drug.map(parent).replace({"CHOLINE FENOFIBRATE": "FENOFIBRATE", "BARDOXOLONE METHYL": "BARDOXOLONE", "CERTOLIZUMAB PEGOL": "CERTOLIZUMAB"})
    r["grade"] = r.parent.map(g).fillna(0).astype(int)
    ext = pd.read_csv(EV / "tie_extension_pool.csv")
    ext["parent"] = ext.drug.map(parent).replace({"BARDOXOLONE METHYL": "BARDOXOLONE", "CERTOLIZUMAB PEGOL": "CERTOLIZUMAB"})
    ext["grade"] = ext.parent.map(g).fillna(0).astype(int)
    pool = pd.concat([r[r.method.isin(METHODS) & (r["rank"] <= 20)][["Drug_ID", "drug", "parent", "grade"]], ext[["Drug_ID", "drug", "parent", "grade"]]]).drop_duplicates("Drug_ID")
    ideal = sorted(pool.grade, reverse=True)
    rows = []
    rng = np.random.default_rng(0)
    for m in METHODS:
        t0 = r[r.method == m].copy()
        for k in (10, 20):
            acc = []
            for _ in range(500):                      # ties at the cut are broken at random, results averaged
                t0["_r"] = rng.random(len(t0))
                top = t0.sort_values(["score", "_r"], ascending=[False, True]).head(k)
                acc.append([(top.grade >= 1).mean(), (top.grade >= 2).mean(), (top.grade == 3).sum(), (top.grade == 2).sum(),
                            (top.grade == 1).sum(), (top.grade == 0).sum(), top.parent.isin(HARM).sum(),
                            dcg(top.grade) / dcg(ideal[:k]) if dcg(ideal[:k]) > 0 else np.nan, top.parent.nunique()])
            a = np.nanmean(acc, axis=0)
            sd = np.nanstd([x[0] for x in acc])
            if m == "SP":                              # 233 of 299 drugs share the top score: the top K is arbitrary and mostly unjudged
                a = np.full(len(a), np.nan); sd = np.nan
            rows.append({"method": m, "K": k, "support_any": a[0], "support_sd_over_ties": sd, "support_animal_or_human": a[1],
                         "human": a[2], "animal": a[3], "in_vitro": a[4], "unknown": a[5], "harm_signal": a[6], "ndcg": a[7],
                         "distinct_parents": a[8]})
    out = pd.DataFrame(rows)
    out.to_csv(MAIN / "table_A_evidence.csv", index=False)
    print(out.round(3).to_string(index=False))
    print("pool drugs", len(pool), "pool parents", pool.parent.nunique(), "grades in pool", pool.drop_duplicates("parent").grade.value_counts().to_dict())
    pool[["Drug_ID", "drug", "parent", "grade"]].sort_values("drug").to_csv(EV / "judged_pool.csv", index=False)


def independent():
    """Evidence support of the top-K lists obtained with the independent-cohort signature (same grades, ties at random)."""
    g = pd.read_csv(EV / "evidence_grades.csv").set_index("parent").grade.to_dict()
    alias = {"CHOLINE FENOFIBRATE": "FENOFIBRATE", "BARDOXOLONE METHYL": "BARDOXOLONE", "CERTOLIZUMAB PEGOL": "CERTOLIZUMAB"}
    rows, rng = [], np.random.default_rng(1)
    pool_g = []
    for m in METHODS:
        f = ROOT / f"results/Evaluation/independent/ranked_full_{m}.csv"
        t = pd.read_csv(f)
        t["parent"] = t.drug.map(parent).replace(alias)
        t["grade"] = t.parent.map(g)
        pool_g += list(t.head(20).grade.dropna())
    for m in METHODS:
        if m == "SP":
            continue
        t = pd.read_csv(ROOT / f"results/Evaluation/independent/ranked_full_{m}.csv")
        t["parent"] = t.drug.map(parent).replace(alias)
        t["grade"] = t.parent.map(g)
        for k in (10, 20):
            sc = t.score.iloc[min(k, len(t)) - 1]
            tie_unjudged = int(t[t.score >= sc].grade.isna().sum())
            acc = []
            for _ in range(300):
                t["_r"] = rng.random(len(t))
                top = t.sort_values(["score", "_r"], ascending=[False, True]).head(k)
                gr = top.grade.fillna(0)
                acc.append([(gr >= 1).mean(), (gr >= 2).mean(), top.parent.isin(HARM).sum(), top.grade.isna().sum()])
            a = np.mean(acc, axis=0)
            rows.append({"method": m, "K": k, "n_ranked": len(t), "support_any": a[0], "support_animal_or_human": a[1],
                         "harm_signal": a[2]})
    out = pd.DataFrame(rows)
    out.to_csv(MAIN / "table_independent_evidence.csv", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
    independent()
