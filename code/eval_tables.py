"""Collect the evaluation results into the tables of the paper (results/Evaluation/tables/*.csv and a printed summary).

Needs the outputs of eval_run.py, eval_robust.py, eval_cost.py, eval_independent.py and eval_metrics.py.

    python code/eval_tables.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from evalcore import ROOT, Data, adjacency, pathway_links  # noqa: E402
from eval_run import all_methods, jac, ranked, topk  # noqa: E402
from eval_evidence_fetch import parent  # noqa: E402

EV = ROOT / "results/Evaluation"
OUT = EV / "tables"
OUT.mkdir(exist_ok=True)
ALIAS = {"CHOLINE FENOFIBRATE": "FENOFIBRATE", "BARDOXOLONE METHYL": "BARDOXOLONE", "CERTOLIZUMAB PEGOL": "CERTOLIZUMAB"}
HARM = {"DAUNORUBICIN", "DOXORUBICIN", "FLOXURIDINE", "DOXIFLURIDINE", "CERTOLIZUMAB", "GOLIMUMAB", "TEGAFUR", "TOPIRAMATE"}


def dcg(rel):
    rel = np.asarray(rel, dtype=float)
    return float(((2 ** rel - 1) / np.log2(np.arange(len(rel)) + 2)).sum())


def evidence(D, s, grades, ideal, ks=(10, 20), draws=500, seed=0):
    t0 = ranked(s).assign(parent=lambda t: t.Drug_ID.map(D.name).map(parent).replace(ALIAS))
    t0["grade"] = t0.parent.map(grades).fillna(0)
    rng = np.random.default_rng(seed)
    out = {}
    for k in ks:
        acc = []
        for _ in range(draws):
            t0["_r"] = rng.random(len(t0))
            top = t0.sort_values(["score", "_r"], ascending=[False, True]).head(k)
            acc.append([(top.grade >= 1).mean(), dcg(top.grade) / dcg(ideal[:k]), (top.grade == 3).sum(), (top.grade == 2).sum(),
                        (top.grade == 1).sum(), top.parent.isin(HARM).sum()])
        a = np.mean(acc, axis=0)
        out[k] = dict(zip(["support", "ndcg", "human", "animal", "in_vitro", "harm"], a))
    return out


def main():
    D = Data()
    A = adjacency(D)
    linked = pathway_links(D, D.seeds)
    sc, cost, R = all_methods(D, A, linked, D.seeds.i.to_numpy())
    grades = pd.read_csv(EV / "evidence/evidence_grades.csv").set_index("parent").grade.to_dict()
    pool = pd.read_csv(EV / "evidence/judged_pool.csv")
    ideal = sorted(pool.grade, reverse=True)
    names = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed", "Decidable-only", "Consistent-any", "Direction-aware DWPC"]
    rows = []
    for m in names:
        ev = evidence(D, sc[m], grades, ideal)
        for k in (10, 20):
            rows.append({"method": m, "K": k, **ev[k], "ranked": int(((sc[m] > 0) if sc[m].min() > -1e8 else (sc[m] > -1e8)).sum())})
    A_tab = pd.DataFrame(rows)
    A_tab.loc[A_tab.method == "SP", ["support", "ndcg", "human", "animal", "in_vitro", "harm"]] = np.nan
    A_tab.to_csv(OUT / "evidence_all_configs.csv", index=False)
    print(A_tab.round(3).to_string(index=False))

    rob = pd.read_csv(EV / "robust/robustness_long.csv")
    st = rob[rob.kind.isin(["seeds", "edges", "noop"]) & (rob.K == 10)].groupby(["kind", "level", "method"]).jaccard.agg(["mean", "std"]).reset_index()
    st.to_csv(OUT / "stability_K10.csv", index=False)
    st20 = rob[rob.kind.isin(["seeds", "edges", "noop"]) & (rob.K == 20)].groupby(["kind", "level", "method"]).jaccard.agg(["mean", "std"]).reset_index()
    st20.to_csv(OUT / "stability_K20.csv", index=False)
    hub = rob[rob.kind.isin(["hub", "param"])].assign(sensitivity=lambda d: 1 - d.jaccard)
    hub.to_csv(OUT / "hub_param_sensitivity.csv", index=False)
    piv = lambda k, lv: st[(st.kind == k) & (st.level == lv)].set_index("method")["mean"]  # noqa: E731
    print("noop ceiling K10", piv("noop", 0.0).round(2).to_dict())
    print("seeds 10% K10", piv("seeds", 0.10).round(2).to_dict())
    print("edges 10% K10", piv("edges", 0.10).round(2).to_dict())
    print(hub[(hub.K == 10)].pivot_table(index=["kind", "level"], columns="method", values="sensitivity").round(2).to_string())


if __name__ == "__main__":
    main()
