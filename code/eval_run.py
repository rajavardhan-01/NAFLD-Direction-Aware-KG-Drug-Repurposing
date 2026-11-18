"""Main comparison: proposed method vs six baselines, ablation family, overlap and cost.

Reads results/Evaluation/inputs (code/eval_export_graph.py). Writes results/Evaluation/main/*.csv.

    python code/eval_run.py
"""
import json
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from evalcore import (ROOT, Data, Routes, adjacency, hetesim_scores, nearest_seed_excl_self, np_scores, path_scores,  # noqa: E402
                      pathway_links, rwr_scores, sp_scores, start_links)

OUT = ROOT / "results/Evaluation/main"
OUT.mkdir(parents=True, exist_ok=True)
KS = (10, 20)


def valid_mask(s):
    return s > -1e8 if s.min() <= -1e8 else s > 0


def ranked(s, rng=None):
    """Drugs with support, best first. Ties are broken by drug id (deterministic) or, with `rng`, at random."""
    m = valid_mask(s)
    t = s[m].reset_index()
    t.columns = ["Drug_ID", "score"]
    if rng is None:
        return t.sort_values(["score", "Drug_ID"], ascending=[False, True]).reset_index(drop=True)
    t["_r"] = rng.random(len(t))
    return t.sort_values(["score", "_r"], ascending=[False, True]).drop(columns="_r").reset_index(drop=True)


def topk(s, k, rng=None):
    return list(ranked(s, rng).Drug_ID[:k])


def tie_size(s, k):
    """Number of drugs that share the score of the K-th ranked drug (1 = no tie at the cut)."""
    t = ranked(s)
    if len(t) < k:
        return 0
    return int((t.score == t.score.iloc[k - 1]).sum())


def jac_ties(sa, sb, k, draws=50, seed=0, fam=None):
    """Mean Jaccard of the top-K of two score vectors over random tie-breaking (independent draws for each)."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(draws):
        a, b = topk(sa, k, rng), topk(sb, k, rng)
        if fam is not None:
            a, b = {fam.get(d, (d,)) for d in a}, {fam.get(d, (d,)) for d in b}
        out.append(jac(a, b))
    return float(np.nanmean(out))


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else np.nan


def all_methods(D, A, linked, seeds_i, alpha=0.4, r=0.5):
    """Score every method once on the full graph. Returns dict name -> Series (and cost per method)."""
    scores, cost = {}, {}

    def timed(name, fn):
        tracemalloc.start()
        t0 = time.perf_counter()
        v = fn()
        dt = time.perf_counter() - t0
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        scores[name], cost[name] = v, {"runtime_s": dt, "peak_MB": peak / 1e6}

    nd = {}
    def _nd():
        nd["v"] = nearest_seed_excl_self(A, seeds_i)
        return nd["v"]
    timed("_dist", _nd)
    timed("SP", lambda: sp_scores(D, nd["v"]))
    timed("NP", lambda: np_scores(D, nd["v"], A))
    for rr in (0.3, 0.5, 0.7):
        timed(f"RWR(r={rr})", lambda rr=rr: rwr_scores(D, A, seeds_i, r=rr))
    box = {}
    def _routes():
        box["R"] = Routes(D, A, seeds_i, linked)
        return None
    timed("_routes", _routes)
    R = box["R"]
    allv = start_links(D, seeds_i, "none")
    prop = start_links(D, seeds_i, "proposed")
    timed("MPC", lambda: path_scores(D, R, allv))
    for a in (0.2, 0.4, 0.6):
        timed(f"DWPC(a={a})", lambda a=a: path_scores(D, R, allv, alpha=a))
    timed("HeteSim", lambda: hetesim_scores(D, R, allv))
    timed("Proposed", lambda: path_scores(D, R, prop))
    # ablation family
    timed("Decidable-only", lambda: path_scores(D, R, start_links(D, seeds_i, "decidable")))
    timed("Consistent-any", lambda: path_scores(D, R, start_links(D, seeds_i, "consistent_any")))
    timed("Direction-aware DWPC", lambda: path_scores(D, R, prop, alpha=alpha))
    # shared setup cost is added to the methods that need it
    cost["SP"]["runtime_s"] += cost["_dist"]["runtime_s"]; cost["NP"]["runtime_s"] += cost["_dist"]["runtime_s"]
    for k in ("MPC", "DWPC(a=0.4)", "HeteSim", "Proposed", "Decidable-only", "Consistent-any", "Direction-aware DWPC"):
        cost[k]["runtime_s"] += cost["_routes"]["runtime_s"]
    scores["DWPC"], cost["DWPC"] = scores["DWPC(a=0.4)"], cost["DWPC(a=0.4)"]
    scores["RWR"], cost["RWR"] = scores["RWR(r=0.5)"], cost["RWR(r=0.5)"]
    return scores, cost, R


def main():
    D = Data()
    A = adjacency(D)
    seeds_i = D.seeds.i.to_numpy()
    linked = pathway_links(D, D.seeds)
    assert linked == D.nafld_pathways
    scores, cost, R = all_methods(D, A, linked, seeds_i)
    main_m = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]
    # --- rankings and overlap
    rk = {m: ranked(scores[m]) for m in scores if not m.startswith("_")}
    pd.concat({m: r.assign(rank=np.arange(1, len(r) + 1)) for m, r in rk.items()}).reset_index(level=0).rename(
        columns={"level_0": "method"}).assign(drug=lambda d: d.Drug_ID.map(D.name)).to_csv(OUT / "all_rankings.csv", index=False)
    rows = []
    for k in KS:
        for a in main_m:
            for b in main_m:
                rows.append({"K": k, "a": a, "b": b, "jaccard": jac_ties(scores[a], scores[b], k)})
    pd.DataFrame(rows).to_csv(OUT / "overlap_jaccard.csv", index=False)
    pd.DataFrame([{"method": m, "K": k, "tied_at_cut": tie_size(scores[m], k), "ranked": int(valid_mask(scores[m]).sum())}
                  for m in main_m for k in KS]).to_csv(OUT / "ties_at_cut.csv", index=False)
    # --- coverage and cost
    cov = []
    for m in main_m + ["Decidable-only", "Consistent-any", "Direction-aware DWPC", "RWR(r=0.3)", "RWR(r=0.7)", "DWPC(a=0.2)", "DWPC(a=0.6)"]:
        c = cost[m]
        cov.append({"method": m, "drugs_with_support": int(valid_mask(scores[m]).sum()), "universe": len(D.universe),
                    "runtime_s": round(c["runtime_s"], 3), "peak_MB": round(c["peak_MB"], 1)})
    cov = pd.DataFrame(cov)
    cov.to_csv(OUT / "coverage_cost.csv", index=False)
    print(cov.to_string(index=False))
    # --- ablation attrition
    att = []
    for mode in ("none", "decidable", "consistent_any", "proposed"):
        use = start_links(D, seeds_i, mode)
        s = path_scores(D, R, use)
        att.append({"filter": mode, "drugs_retained": use.Drug_ID.nunique(), "start_links": len(use),
                    "drugs_with_BB_routes": int((s > 0).sum()), "BB_routes": int(s.sum())})
    pd.DataFrame(att).to_csv(OUT / "ablation_attrition.csv", index=False)
    print(pd.DataFrame(att).to_string(index=False))
    # --- ablation overlap with proposed
    ab = ["MPC", "Decidable-only", "Consistent-any", "Proposed", "DWPC", "Direction-aware DWPC"]
    rows = [{"K": k, "config": a, "jaccard_with_proposed": jac_ties(scores[a], scores["Proposed"], k)}
            for k in KS for a in ab]
    pd.DataFrame(rows).to_csv(OUT / "ablation_overlap.csv", index=False)
    print(pd.DataFrame(rows).pivot(index="config", columns="K", values="jaccard_with_proposed").round(3))
    # sensitivity of RWR / DWPC parameters to the default (self-overlap)
    sens = []
    for k in KS:
        for a, b in [("RWR(r=0.5)", "RWR(r=0.3)"), ("RWR(r=0.5)", "RWR(r=0.7)"), ("DWPC(a=0.4)", "DWPC(a=0.2)"), ("DWPC(a=0.4)", "DWPC(a=0.6)")]:
            sens.append({"K": k, "default": a, "variant": b, "jaccard": jac_ties(scores[a], scores[b], k)})
    pd.DataFrame(sens).to_csv(OUT / "parameter_sensitivity.csv", index=False)
    print(pd.DataFrame(sens).to_string(index=False))
    info = {"universe": len(D.universe), "np_zero_variance": int(scores["NP"].attrs.get("zero_variance", -1)),
            "np_no_reachable_target": int(scores["NP"].attrs.get("no_reachable_target", -1)),
            "top20_union": len({d for m in main_m for d in topk(scores[m], 20)}),
            "top10_union": len({d for m in main_m for d in topk(scores[m], 10)})}
    (OUT / "info.json").write_text(json.dumps(info, indent=2))
    print(info)


if __name__ == "__main__":
    main()
