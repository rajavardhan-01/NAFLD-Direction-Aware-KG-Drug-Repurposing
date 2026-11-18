"""Robustness study: seed removal, edge removal, hub removal and parameter sensitivity for all seven methods.

For every random perturbation level the SAME perturbed graph is given to all methods (shared perturbation) and
the seed-derived pathway selection is recomputed. Stability = Jaccard of the top-K drugs against the unperturbed
ranking (also at target-family level). Writes results/Evaluation/robust/*.csv.

    python code/eval_robust.py [reps]          # default 100 repetitions per level
"""
import copy
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from evalcore import (HUB_MAX, ROOT, Data, Routes, adjacency, hetesim_scores, nearest_seed_excl_self, np_scores,  # noqa: E402
                      path_scores, pathway_links, rwr_scores, sp_scores, start_links)
from eval_run import jac, topk  # noqa: E402
DRAWS_DET = 20

OUT = ROOT / "results/Evaluation/robust"
KS = (10, 20)
METHODS = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed", "Decidable-only", "Direction-aware DWPC"]
_D = None
_BASE = None


def score_main(D, A, seeds_i, linked, memb=None, hub=HUB_MAX, links_override=None):
    if links_override is not None:
        D = copy.copy(D)
        D.links = links_override
    nd = nearest_seed_excl_self(A, seeds_i)
    R = Routes(D, A, seeds_i, linked, memb=memb, hub=hub)
    allv = start_links(D, seeds_i, "none")
    prop = start_links(D, seeds_i, "proposed")
    return {"SP": sp_scores(D, nd), "NP": np_scores(D, nd, A), "RWR": rwr_scores(D, A, seeds_i, r=0.5),
            "MPC": path_scores(D, R, allv), "DWPC": path_scores(D, R, allv, alpha=0.4),
            "HeteSim": hetesim_scores(D, R, allv), "Proposed": path_scores(D, R, prop),
            "Decidable-only": path_scores(D, R, start_links(D, seeds_i, "decidable")),
            "Direction-aware DWPC": path_scores(D, R, prop, alpha=0.4)}


def families(D):
    """Target family of a drug = the set of its seed-linked target genes (drugs with the same set act through one mechanism)."""
    return D.links.groupby("Drug_ID").Gene_ID.agg(lambda s: tuple(sorted(set(s)))).to_dict()


def fam_set(fam, drugs):
    return {fam.get(d, (d,)) for d in drugs}


def init():
    global _D, _BASE
    _D = Data()
    A = adjacency(_D)
    linked = pathway_links(_D, _D.seeds)
    seeds_i = _D.seeds.i.to_numpy()
    sc = score_main(_D, A, seeds_i, linked)
    _BASE = {"A": A, "linked": linked, "seeds_i": seeds_i, "scores": sc, "fam": families(_D),
             "top": {m: {k: topk(sc[m], k) for k in KS} for m in METHODS}}


def compare(sc, tag, draws=1, seed=0):
    """Top-K overlap with the unperturbed ranking; tied scores are broken at random (independent draws)."""
    rows = []
    D, B = _D, _BASE
    rng = np.random.default_rng(seed)
    for m in METHODS:
        s = sc[m]
        mask = s > -1e8 if s.min() <= -1e8 else s > 0
        for k in KS:
            j, f = [], []
            for _ in range(draws):
                t = topk(s, k, rng)
                o = topk(B["scores"][m], k, rng)
                j.append(jac(o, t))
                f.append(jac(fam_set(B["fam"], o), fam_set(B["fam"], t)))
            rows.append({**tag, "method": m, "K": k, "jaccard": float(np.nanmean(j)), "family_jaccard": float(np.nanmean(f)),
                         "n_support": int(mask.sum())})
    return rows


def perturb_task(args):
    kind, level, rep = args
    D, B = _D, _BASE
    rng = np.random.default_rng(1000 * rep + int(level * 100) + (0 if kind == "seeds" else 7))
    seeds, A, memb, links = D.seeds, B["A"], D.memb, D.links
    Dm = D
    if kind == "seeds":
        keep = []
        for _, g in seeds.groupby("direction"):
            n_rm = int(round(level * len(g)))
            rm = rng.choice(g.index.to_numpy(), size=n_rm, replace=False)
            keep += [i for i in g.index if i not in set(rm)]
        seeds = seeds.loc[sorted(keep)]
    elif kind == "noop":
        pass
    elif kind == "edges":
        e = D.edges
        e2 = e[rng.random(len(e)) >= level]
        A = adjacency(D, edges=e2)
        dt = D.drug_targets_all
        gone = rng.random(len(dt)) < level
        dt2 = dt[~gone]
        key_gone = set(zip(dt.drug_id[gone], dt.i[gone]))
        links = D.links[[(a, b) not in key_gone for a, b in zip(D.links.Drug_ID, D.links.i)]]
        memb = D.memb[rng.random(len(D.memb)) >= level].reset_index(drop=True)
        Dm = copy.copy(D)
        Dm.drug_targets_all, Dm.links, Dm.memb = dt2, links, memb
        Dm.gene_indeg = np.bincount(dt2.i, minlength=D.N)
        Dm.drug_outdeg = dt2.groupby("drug_id").size()
    seeds_i = seeds.i.to_numpy()
    linked = pathway_links(Dm, seeds, memb if kind != "seeds" else None) if kind in ("seeds", "edges") else B["linked"]
    sc = score_main(Dm, A, seeds_i, linked, memb=memb if kind == "edges" else None)
    return compare(sc, {"kind": kind, "level": level, "rep": rep}, seed=rep + 7919 * int(level * 100))


def hub_task(args):
    name, pct, genes = args
    D, B = _D, _BASE
    A = adjacency(D, removed_nodes=np.array(genes, dtype=int))
    sc = score_main(D, A, B["seeds_i"], B["linked"])
    return compare(sc, {"kind": "hub", "level": name, "rep": 0}, draws=DRAWS_DET)


def param_task(args):
    name, kw = args
    D, B = _D, _BASE
    kw = dict(kw)
    links_override = None
    if kw.pop("cleavage", False):
        l = D.links.copy()
        raw = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_links_direction_check.csv", usecols=["Drug_ID", "Gene_ID", "Action_raw"])
        hyd = set(zip(raw.Drug_ID[raw.Action_raw == "HYDROLYTIC_ENZYME"], raw.Gene_ID[raw.Action_raw == "HYDROLYTIC_ENZYME"]))
        l.loc[[(a, b) in hyd for a, b in zip(l.Drug_ID, l.Gene_ID)], "dcheck"] = "undecided"
        links_override = l
    sc = score_main(D, B["A"], B["seeds_i"], B["linked"], links_override=links_override, **kw)
    return compare(sc, {"kind": "param", "level": name, "rep": 0}, draws=DRAWS_DET)


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    init()
    A = _BASE["A"]
    deg = np.asarray(A.sum(1)).ravel()
    ppi_nodes = np.where(deg > 0)[0]
    order = ppi_nodes[np.argsort(-deg[ppi_nodes], kind="stable")]
    ppara = _D.idx["ENSG00000186951"]
    hub_jobs = [("top1%", 1, order[: int(round(0.01 * len(order)))].tolist()), ("top5%", 5, order[: int(round(0.05 * len(order)))].tolist()),
                ("PPARA", 0, [ppara])]
    param_jobs = [(f"hub<={h}", {"hub": h}) for h in (100, 250, 1000, 2000, 10 ** 9)] + [("cleavage->undecided", {"cleavage": True})]
    jobs = [("seeds", f, r) for f in (0.05, 0.10, 0.15, 0.20) for r in range(reps)] + \
           [("edges", f, r) for f in (0.05, 0.10, 0.20) for r in range(reps)] +            [("noop", 0.0, r) for r in range(reps)]
    print("PPI nodes", len(order), "hub sizes", [len(h[2]) for h in hub_jobs], "jobs", len(jobs), flush=True)
    rows = []
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 8, initializer=init) as pool:
        for r in pool.imap_unordered(hub_task, hub_jobs):
            rows += r
        for r in pool.imap_unordered(param_task, param_jobs):
            rows += r
        done = 0
        for r in pool.imap_unordered(perturb_task, jobs, chunksize=2):
            rows += r
            done += 1
            if done % 50 == 0:
                print(f"{done}/{len(jobs)} {time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "robustness_long.csv", index=False)
    summ = (df[df.kind.isin(["seeds", "edges", "noop"])].groupby(["kind", "level", "method", "K"])
            .agg(jaccard_mean=("jaccard", "mean"), jaccard_sd=("jaccard", "std"), jaccard_lo=("jaccard", lambda x: x.quantile(0.025)),
                 jaccard_hi=("jaccard", lambda x: x.quantile(0.975)), family_mean=("family_jaccard", "mean"),
                 support_mean=("n_support", "mean"), reps=("rep", "count")).reset_index())
    summ.to_csv(OUT / "robustness_summary.csv", index=False)
    print(summ[summ.K == 10].pivot_table(index=["kind", "level"], columns="method", values="jaccard_mean").round(2))
    print(df[df.kind.isin(["hub", "param"]) & (df.K == 10)].pivot_table(index=["kind", "level"], columns="method", values="jaccard").round(2))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
