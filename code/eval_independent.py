"""Independent-cohort evaluation: a new NAFLD liver signature from GEO GSE126848 (Suppli et al. 2019, liver biopsies),
run through the unchanged proposed method and the six baselines.

Cohort: healthy normal-weight (n=14) and obese (n=12) controls, simple steatosis NAFL (n=15), NASH (n=16); RNA-seq counts.
Signature (thresholds fixed before the analysis): log2(DESeq-normalised count + 1) ~ group + sex, NAFLD (NAFL+NASH) versus
all controls; kept if BH-FDR < 0.01 and |log2 fold change| >= 1. Stage annotation: 2_NASH if the gene is also significant
in NASH versus controls (same sign), otherwise 1_steatosis. Direction = sign of the fold change.
The pathway links are rebuilt with the same enrichment + CTD rule; the candidate universe (299 drugs) is held fixed and
its drug-gene links to the new seeds are read from the graph with the same action-to-direction rule.

    python code/eval_independent.py
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist

sys.path.insert(0, str(Path(__file__).parent))
from evalcore import ROOT, Data, Routes, adjacency, pathway_links, start_links  # noqa: E402
from eval_run import jac, jac_ties, ranked, topk  # noqa: E402
from eval_robust import families, fam_set, score_main  # noqa: E402

GEO = ROOT / "data/geo"
OUT = ROOT / "results/Evaluation/independent"
OUT.mkdir(parents=True, exist_ok=True)
FDR, LFC = 0.01, 1.0
DOWN = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE_MODULATOR", "NEGATIVE_ALLOSTERIC_MODULATOR",
        "ALLOSTERIC_ANTAGONIST", "INVERSE_AGONIST", "DEGRADER"}
UP = {"AGONIST", "PARTIAL_AGONIST", "ACTIVATOR", "POSITIVE_MODULATOR", "POSITIVE_ALLOSTERIC_MODULATOR", "OPENER", "STABILISER"}
MIXED = {"MODULATOR", "BINDING_AGENT", "SUBSTRATE", "RELEASING_AGENT"}
PUSH = {"down": "down", "up": "up", "conflicting": "mixed", "unknown": "unknown"}


def adir(a):
    a = str(a).upper()
    return "down" if a in DOWN else "up" if a in UP else "mixed" if a in MIXED else "unknown"


def blank(x):
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == ""


def link_direction(row):
    """Action direction of a link, the rule of code/task2_drug_gene_edges.py."""
    added = row.layer == "added curated drug-target links"
    if added:
        return PUSH.get(str(row.ext_push), "unknown")
    d = adir(row.okg_action)
    ass = row.relation_assertions
    if isinstance(ass, str) and ass and {"down", "up"} <= {adir(x.split(":", 1)[-1]) for x in ass.split(";")}:
        d = "mixed"
    if d == "unknown" and row.external_action_direction in ("down", "up"):
        d = row.external_action_direction
    return d


def read_counts():
    c = pd.read_csv(GEO / "GSE126848_Gene_counts_raw.txt.gz", sep="\t", index_col=0)
    t = gzip.open(GEO / "GSE126848_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace").read().split("\n")
    get = lambda key: [x.strip('"') for x in next(l for l in t if l.startswith(key)).split("\t")[1:]]  # noqa: E731
    desc = get("!Sample_description")
    chars = [[x.strip('"') for x in l.split("\t")[1:]] for l in t if l.startswith("!Sample_characteristics_ch1")]
    meta = pd.DataFrame({"col": desc, "sex": [x.split(": ")[1] for x in chars[0]], "group": [x.split(": ")[1] for x in chars[2]]})
    c.columns = c.columns.astype(str)
    meta["col"] = meta.col.astype(str).str.zfill(4)
    assert set(meta.col) <= set(c.columns), (meta.col[:5].tolist(), c.columns[:5].tolist())
    return c[meta.col], meta


def deseq_norm(c):
    lc = np.log(c.replace(0, np.nan))
    ref = lc.mean(axis=1)
    ok = ref.notna() & np.isfinite(ref)
    sf = np.exp((lc[ok].sub(ref[ok], axis=0)).median(axis=0))
    return c / sf


def fit(y, X, k):
    """OLS of every row of y on X; returns coefficient k and its two-sided p-value."""
    XtXi = np.linalg.inv(X.T @ X)
    beta = y @ X @ XtXi
    res = y - beta @ X.T
    df = X.shape[0] - X.shape[1]
    s2 = (res ** 2).sum(1) / df
    se = np.sqrt(s2 * XtXi[k, k])
    tt = beta[:, k] / se
    return beta[:, k], 2 * tdist.sf(np.abs(tt), df)


def signature(c, meta):
    c = c[(c >= 10).sum(axis=1) >= 12]
    y = np.log2(deseq_norm(c) + 1).to_numpy()
    sex = (meta.sex == "Male").astype(float).to_numpy()
    nafld = meta.group.isin(["NAFLD", "NASH"]).astype(float).to_numpy()
    nash = (meta.group == "NASH").astype(float).to_numpy()
    ctrl = meta.group.isin(["healthy", "obese"]).to_numpy()
    X = np.c_[np.ones(len(meta)), nafld, sex]
    lfc, p = fit(y, X, 1)
    sel = ctrl | (meta.group == "NASH").to_numpy()
    Xn = np.c_[np.ones(sel.sum()), nash[sel], sex[sel]]
    lfc_n, p_n = fit(y[:, sel], Xn, 1)
    from evalcore import bh_qvalues
    q, qn = bh_qvalues(p), bh_qvalues(p_n)
    d = pd.DataFrame({"gene_id": c.index, "log2FC": lfc, "p": p, "q": q, "log2FC_NASH": lfc_n, "q_NASH": qn})
    d["significant"] = (d.q < FDR) & (d.log2FC.abs() >= LFC)
    d["stage"] = np.where((d.q_NASH < 0.05) & (np.sign(d.log2FC_NASH) == np.sign(d.log2FC)), "2_NASH", "1_steatosis")
    d["direction"] = np.where(d.log2FC > 0, "increasing", "decreasing")
    return d


def rbo(a, b, p=0.9, depth=20):
    """Rank-biased overlap (truncated at `depth`, extrapolated)."""
    a, b = list(a)[:depth], list(b)[:depth]
    depth = min(len(a), len(b))
    if depth == 0:
        return float('nan')
    s, A, B = 0.0, set(), set()
    for d in range(1, depth + 1):
        A.add(a[d - 1]); B.add(b[d - 1])
        s += (p ** (d - 1)) * len(A & B) / d
    x = len(A & B) / depth
    return (1 - p) * s + x * p ** depth


def links_for(D, links_all, seeds):
    L = links_all.copy()
    L = L[L.Drug_ID.isin(D.universe) & L.Gene_ID.isin(set(seeds.ensembl_id))]
    L["Action_direction"] = L.apply(link_direction, axis=1)
    sd = seeds.set_index("ensembl_id").direction
    gd = L.Gene_ID.map(sd)
    a = L.Action_direction
    L["dcheck"] = np.select([((a == "down") & (gd == "increasing")) | ((a == "up") & (gd == "decreasing")),
                             ((a == "down") & (gd == "decreasing")) | ((a == "up") & (gd == "increasing"))], ["consistent", "inconsistent"], "undecided")
    L["i"] = L.Gene_ID.map(D.idx)
    L["Drug"] = L.Drug_ID.map(D.name)
    return L[["Drug", "Drug_ID", "Gene_ID", "i", "dcheck"]].drop_duplicates().reset_index(drop=True)


def main():
    D = Data()
    A = adjacency(D)
    la = pd.read_csv(ROOT / "results/Evaluation/inputs/candidate_links_all.csv")
    # check: the direction rule reproduces the stored direction check on the original seeds
    chk = links_for(D, la, D.seeds.assign(direction=D.seeds.direction))
    ref = D.links[D.links.Drug_ID.isin(D.universe)].set_index(["Drug_ID", "Gene_ID"]).dcheck
    mm = chk.set_index(["Drug_ID", "Gene_ID"]).dcheck.reindex(ref.index)
    print("rule check on original seeds: links", len(ref), "mismatches", int((mm != ref).sum()))
    c, meta = read_counts()
    print(meta.group.value_counts().to_dict(), c.shape)
    sig = signature(c, meta)
    sig.to_csv(OUT / "differential_expression.csv", index=False)
    s = sig[sig.significant & sig.gene_id.isin(D.idx)].copy()
    n_all = int(sig.significant.sum())
    s["ensembl_id"] = s.gene_id
    s["i"] = s.gene_id.map(D.idx)
    s["main_stage"] = s.stage
    orig = D.seeds
    print("DE tested", len(sig), "significant", n_all, "mapped to graph", len(s), "up", int((s.direction == "increasing").sum()))
    ov = set(s.ensembl_id) & set(orig.ensembl_id)
    conc = np.mean([(s.set_index("ensembl_id").direction[g] == orig.set_index("ensembl_id").direction[g]) for g in ov
                    if orig.set_index("ensembl_id").direction[g] in ("increasing", "decreasing")])
    print("overlap with original seeds", len(ov), "direction concordance", round(conc, 3))
    variants = {"full signature": s,
                "size-matched (top 280 by FDR)": s.sort_values("q").head(280)}
    base_linked = pathway_links(D, orig)
    scores0 = score_main(D, A, orig.i.to_numpy(), base_linked)
    fam = families(D)
    rows, info = [], []
    for name, seeds in variants.items():
        seeds = seeds.reset_index(drop=True)
        linked = pathway_links(D, seeds)
        L = links_for(D, la, seeds)
        D2 = __import__("copy").copy(D)
        D2.links = L
        sc = score_main(D2, A, seeds.i.to_numpy(), linked)
        R = Routes(D2, A, seeds.i.to_numpy(), linked)
        info.append({"signature": name, "seeds": len(seeds), "increasing": int((seeds.direction == "increasing").sum()),
                     "decreasing": int((seeds.direction == "decreasing").sum()), "stage_NASH": int((seeds.main_stage == "2_NASH").sum()),
                     "seed_overlap_with_original": len(set(seeds.ensembl_id) & set(orig.ensembl_id)),
                     "pathways_linked": len(linked), "pathways_shared_with_original": len(linked & base_linked),
                     "eligible_pathways": len(R.elig),
                     "candidate_links_to_seeds": len(L), "consistent": int((L.dcheck == "consistent").sum()),
                     "inconsistent": int((L.dcheck == "inconsistent").sum()), "undecided": int((L.dcheck == "undecided").sum()),
                     "drugs_direction_consistent": int(start_links(D2, seeds.i.to_numpy(), "proposed").Drug_ID.nunique())})
        for m in ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]:
            for k in (10, 20):
                rng = np.random.default_rng(k)
                rb = np.nanmean([rbo(topk(scores0[m], 20, rng), topk(sc[m], 20, rng)) for _ in range(50)]) if k == 20 else np.nan
                rows.append({"signature": name, "method": m, "K": k, "jaccard": jac_ties(scores0[m], sc[m], k),
                             "family_jaccard": jac_ties(scores0[m], sc[m], k, fam=fam), "rbo": rb,
                             "support_new": int((sc[m] > -1e8).sum() if sc[m].min() <= -1e8 else (sc[m] > 0).sum())})
            tt = ranked(sc[m])
            ext_ids = tt[tt.score >= tt.score.iloc[min(19, len(tt) - 1)]].Drug_ID.tolist() if len(tt) else []
            tt.assign(drug=lambda d: d.Drug_ID.map(D.name)).to_csv(OUT / f"ranked_{name.split()[0]}_{m}.csv", index=False)
            pd.DataFrame({"Drug_ID": ext_ids}).assign(drug=lambda d: d.Drug_ID.map(D.name), method=m, signature=name).to_csv(
                OUT / f"top20ext_{name.split()[0]}_{m}.csv", index=False)
            pd.DataFrame({"Drug_ID": topk(sc[m], 20)}).assign(drug=lambda d: d.Drug_ID.map(D.name), method=m, signature=name).to_csv(
                OUT / f"top20_{name.split()[0]}_{m}.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "cross_signature_stability.csv", index=False)
    pd.DataFrame(info).to_csv(OUT / "signature_summary.csv", index=False)
    print(pd.DataFrame(info).T.to_string())
    print(pd.DataFrame(rows).pivot_table(index=["signature", "method"], columns="K", values=["jaccard", "family_jaccard"]).round(2).to_string())
    (OUT / "overlap_note.txt").write_text(f"DE tested {len(sig)}, significant {n_all}, mapped {len(s)}, overlap {len(ov)}, direction concordance {conc:.3f}\n")


if __name__ == "__main__":
    main()
