"""Step 0.7: Pathway -[:IMPLICATED_IN]-> NAFLD (EFO_0003095) from seed-gene enrichment + CTD curation.

Enrichment: one-sided hypergeometric test of seed genes in Reactome (OptimusKG), KEGG legacy and Hallmark
pathways (MSigDB 2024.1). Background N = all genes with >= 1 IN_PATHWAY link in this graph (any database).
Gene groups (seed file, main_stage and direction columns):
  all_histology            main_stage != 0_genetic_risk
  <stage>                  main_stage == 1_steatosis | 2_NASH | 3_fibrosis | 4_HCC
  <stage>_<direction>      the same, direction == increasing | decreasing
Only seed genes in the background count. BH-FDR is applied within each group x database over every
pathway of that database with >= 1 background gene. Significant = q < 0.05.

Link rule: a pathway is linked to NAFLD only if it is significant in >= 1 group AND
  * Reactome / KEGG: CTD (CTD_diseases_pathways, MESH:D065626) links it to NAFLD through >= 2 curated genes;
  * Hallmark: no CTD condition (CTD has no Hallmark pathways, so the condition can never be met).
The Hallmark exception and the lower minimum (MIN_LINKS = 12) were approved by the user on 2026-11-14 after
the strict two-source rule gave 6 links.
CTD matching: Reactome by stable id (REACT:R-HSA-...); KEGG by normalised pathway name.
Stops (writes nothing) if fewer than MIN_LINKS pathways pass. Dry run by default; --apply writes; --undo removes.
Also writes each seed gene's most significant pathway (used by Task 4).

    python code/pathway_nafld_links.py
    python code/pathway_nafld_links.py --apply
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import NAFLD, ROOT, counts, driver, now, read_session, write_log, write_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
CTD = ROOT / "data/ctd/CTD_diseases_pathways.csv.gz"
OUT = ROOT / "results/graph_build/step0_7_pathway_nafld"
SOURCE = "seed enrichment + CTD"
STAGES = ["1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]
MIN_LINKS = 12   # user-approved 2026-11-14 (was 20)


def db_of(node_id, source):
    if str(node_id).startswith("REACT:"):
        return "Reactome"
    return {"KEGG_legacy": "KEGG", "Hallmark": "Hallmark"}.get(source, "other")


def norm(v):
    return re.sub(r"[^A-Z0-9]+", "_", str(v).upper()).strip("_")


def groups(seeds):
    h = seeds[seeds.main_stage != "0_genetic_risk"]
    out = {"all_histology": set(h.ensembl_id)}
    for st in STAGES:
        sub = seeds[seeds.main_stage == st]
        out[st] = set(sub.ensembl_id)
        for d in ["increasing", "decreasing"]:
            out[f"{st}_{d}"] = set(sub[sub.direction == d].ensembl_id)
    return out


def enrichment(memb, seeds):
    bg = set(memb.gene)
    N = len(bg)
    pw = memb.groupby(["pathway", "name", "database"]).gene.agg(set).reset_index()
    rows = []
    for gname, gset in groups(seeds).items():
        gin = gset & bg
        n = len(gin)
        if n == 0:
            continue
        for db, sub in pw.groupby("database"):
            res = []
            for p, nm, gs in zip(sub.pathway, sub.name, sub.gene):
                k_genes = gs & gin
                K, k = len(gs), len(k_genes)
                pval = float(hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
                res.append((gname, db, p, nm, K, n, N, k, pval, (k / n) / (K / N) if k else 0.0,
                            ";".join(sorted(k_genes))))
            r = pd.DataFrame(res, columns=["group", "database", "pathway", "name", "pathway_size", "group_genes_in_bg",
                                           "background", "overlap", "p_value", "fold_enrichment", "overlap_genes"])
            r["q_value"] = multipletests(r.p_value, method="fdr_bh")[1]
            rows.append(r)
    e = pd.concat(rows, ignore_index=True)
    return e, N


def ctd_support(kegg_nodes):
    c = pd.read_csv(CTD, comment="#", header=None, names=["DiseaseName", "DiseaseID", "PathwayName", "PathwayID", "Genes"])
    c = c[c.DiseaseID.eq("MESH:D065626")]
    sup = {}
    kegg_by_name = {norm(n[len("KEGG:KEGG_"):]): n for n in kegg_nodes}
    for (pid, name), g in c.groupby(["PathwayID", "PathwayName"]):
        genes = {v for t in g.Genes.dropna() for v in re.split(r"[|;]", t)} - {""}
        if pid.startswith("REACT:"):
            sup[pid] = (pid, name, len(genes))
        elif pid.startswith("KEGG:"):
            node = kegg_by_name.get(norm(name))
            if node:
                if node in sup:
                    raise SystemExit(f"ambiguous CTD KEGG name mapping: {name}")
                sup[node] = (pid, name, len(genes))
    return sup, c


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--undo", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    drv = driver()
    if args.undo:
        with write_session(drv) as s:
            n = s.run("MATCH (:Pathway)-[r:IMPLICATED_IN {source: $s}]->(:Disease) DELETE r RETURN count(r)", s=SOURCE).single()[0]
        print("removed", n, "links")
        return
    seeds = pd.read_csv(SEEDS)
    with read_session(drv) as s:
        memb = pd.DataFrame(s.run("MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway) RETURN DISTINCT g.node_id AS gene, "
                                  "p.node_id AS pathway, p.name AS name, p.source AS source").data())
        if not s.run("MATCH (d:Disease {node_id: $n}) RETURN count(d)", n=NAFLD).single()[0]:
            raise SystemExit(f"{NAFLD} not in graph")
        existing = s.run("MATCH ()-[r:IMPLICATED_IN]->() RETURN count(r)").single()[0]
    if existing:
        raise SystemExit("IMPLICATED_IN links already exist; undo first")
    memb["database"] = [db_of(p, s_) for p, s_ in zip(memb.pathway, memb.source)]
    e, N = enrichment(memb, seeds)
    e.sort_values(["group", "database", "p_value"]).to_csv(OUT / "enrichment_all.csv", index=False)
    sig = e[e.q_value < 0.05]
    sig.sort_values(["group", "database", "q_value"]).to_csv(OUT / "enrichment_significant.csv", index=False)
    kegg_nodes = sorted(memb.loc[memb.database == "KEGG", "pathway"].unique())
    sup, ctd = ctd_support(kegg_nodes)

    links, table = [], []
    for p, gp in sig.groupby("pathway"):
        best = gp.sort_values("q_value").iloc[0]
        cs = sup.get(p)
        hallmark = best.database == "Hallmark"
        passes = hallmark or (cs is not None and cs[2] >= 2)
        table.append({"pathway": p, "name": best["name"], "database": best.database, "best_q": best.q_value,
                      "best_group": best.group, "groups": ";".join(sorted(gp.group)), "in_ctd_nafld": cs is not None,
                      "ctd_n_genes": cs[2] if cs else 0, "linked": passes,
                      "link_basis": ("enrichment only (Hallmark, no CTD coverage)" if hallmark else "enrichment + CTD") if passes else ""})
        if passes:
            links.append({"pathway": p, "props": {
                "source": SOURCE, "database": best.database, "best_q": float(best.q_value),
                "best_fold": round(float(best.fold_enrichment), 3), "best_group": best.group,
                "groups": ";".join(sorted(gp.group)), "ctd_pathway_id": cs[0] if cs else "",
                "ctd_n_genes": int(cs[2]) if cs else 0,
                "link_basis": "enrichment only (Hallmark, no CTD coverage)" if hallmark else "enrichment + CTD"}})
    t = pd.DataFrame(table)
    t.to_csv(OUT / "significant_pathways_ctd_check.csv", index=False)

    # each seed gene's most significant pathway (any group containing it; significant only)
    x = sig.assign(gene=sig.overlap_genes.str.split(";")).explode("gene")
    top = x.sort_values(["q_value", "p_value"]).drop_duplicates("gene")
    gene_top = seeds[["gene", "ensembl_id"]].merge(
        top[["gene", "pathway", "name", "database", "q_value", "group"]].rename(
            columns={"gene": "ensembl_id", "name": "pathway_name", "group": "pathway_group"}), on="ensembl_id", how="left")
    gene_top.rename(columns={"pathway": "top_pathway_id", "pathway_name": "top_pathway", "q_value": "top_pathway_q",
                             "database": "top_pathway_db"}).to_csv(OUT / "seed_gene_top_pathway.csv", index=False)

    grp_sizes = {k: len(v & set(memb.gene)) for k, v in groups(seeds).items()}
    summary = {"background_genes": N, "group_sizes_in_background": grp_sizes,
               "pathways_tested": memb.groupby("database").pathway.nunique().to_dict(),
               "significant_pathways": t.groupby("database").size().to_dict() if len(t) else {},
               "ctd_nafld_pathways": {"Reactome": int(ctd.PathwayID.str.startswith("REACT:").sum()),
                                      "KEGG": int(ctd.PathwayID.str.startswith("KEGG:").sum()),
                                      "KEGG_matched_to_GMT": sum(k.startswith("KEGG:") for k in sup)},
               "significant_and_in_ctd": int(t.in_ctd_nafld.sum()) if len(t) else 0,
               "linked": len(links),
               "linked_by_database": pd.Series([l["props"]["database"] for l in links]).value_counts().to_dict()}
    import json
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))
    if len(links) < MIN_LINKS:
        raise SystemExit(f"STOP: only {len(links)} pathways pass (< {MIN_LINKS}); rule not changed, nothing written.")
    if not args.apply:
        print("dry run: nothing written")
        return
    with write_session(drv) as s:
        before = counts(s)
        added = now()
        (OUT / "undo_implicated_in.cypher").write_text(
            f"// Undo step 0.7\nMATCH (:Pathway)-[r:IMPLICATED_IN {{source: '{SOURCE}'}}]->(:Disease) DELETE r;\n",
            encoding="utf-8")

        def tx_write(tx):
            n = tx.run("""UNWIND $rows AS row MATCH (p:Pathway {node_id: row.pathway}) MATCH (d:Disease {node_id: $nafld})
                CREATE (p)-[r:IMPLICATED_IN]->(d) SET r = row.props, r.added = $a RETURN count(r) AS n""",
                       rows=links, nafld=NAFLD, a=added).single()["n"]
            if n != len(links):
                raise RuntimeError("count mismatch; rolled back")
        s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    write_log(OUT / "write_log.json", {"step": "0.7 pathway -> NAFLD", "added": added, "links": len(links),
                                       "summary": summary, "before": before, "after": after,
                                       "undo": "undo_implicated_in.cypher"})
    print(f"written: {len(links)} IMPLICATED_IN links")


if __name__ == "__main__":
    main()
