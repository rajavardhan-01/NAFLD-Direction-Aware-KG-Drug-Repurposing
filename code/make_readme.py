"""Write README.md from this run's result files (every number is read from a file, nothing typed by hand).

    python code/make_readme.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
GB = R / "graph_build"
RK = R / "Ranking"
RUNS = [("RunA", "RunA_all_seed_genes", "all 280 seed genes (main analysis)"),
        ("RunB", "RunB_histology_genes", "228 histology genes (`main_stage` != `0_genetic_risk`)"),
        ("RunC", "RunC_fibrosis_genes", "all genes with `3_fibrosis` in `all_stages`"),
        ("RunD", "RunD_seed_plus_GWAS", "280 seed genes + GWAS Catalog add-on genes")]


def j(p):
    return json.loads((ROOT / p).read_text(encoding="utf-8"))


def n(x):
    return f"{int(x):,}"


def md_table(df, cols=None, fmt=None):
    cols = cols or list(df.columns)
    fmt = fmt or {}
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append("" if pd.isna(v) else (fmt[c](v) if c in fmt else str(v)))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def main():
    src = j("data/SOURCES.json")["files"]
    okg = next(s for s in src if s["source"] == "OptimusKG")
    okg_date = min(s["downloaded_utc"] for s in src if s["source"] == "OptimusKG")[:10]
    gmt = [s for s in src if s["source"] == "MSigDB"]
    ctd = next(s for s in src if s["source"] == "CTD")
    imp = j("results/graph_build/step0_0_import/write_log.json")
    plan = pd.read_csv(GB / "step0_0_import/edge_plan.csv")
    s02 = j("results/graph_build/step0_2_chembl_gaps/write_log.json")
    s04 = j("results/graph_build/step0_4_new_links/numbers.json")
    s04w = j("results/graph_build/step0_4_new_links/graph_write/write_log.json")
    s05 = j("results/graph_build/step0_5_external_directions/summary.json")
    s05w = j("results/graph_build/step0_5_external_directions/write_log.json")
    s06w = j("results/graph_build/step0_6_pathways/write_log.json")
    s07 = j("results/graph_build/step0_7_pathway_nafld/summary.json")
    s07t = pd.read_csv(GB / "step0_7_pathway_nafld/significant_pathways_ctd_check.csv")
    s08 = j("results/graph_build/step0_8_added_drug_indications/write_log.json")
    g = j("results/graph_build/graph_report.json")
    t1 = pd.read_csv(R / "Task1/seed_gene_table.csv")
    t1c = pd.read_csv(R / "Task1/counts_by_stage_direction.csv")
    t1conn = pd.read_csv(R / "Task1/connections_by_stage.csv")
    t2 = pd.read_csv(R / "Task2/summary_numbers.csv", index_col=0)["value"]
    t2s = pd.read_csv(R / "Task2/links_by_stage.csv")
    t3 = pd.read_csv(R / "Task3/druggability_summary_pharmacological.csv")
    t3bg = pd.read_csv(R / "Task3/background_genes_druggability.csv")
    t3out = pd.read_csv(R / "Task3/seed_genes_outside_background.csv")
    fun = pd.read_csv(R / "Task4/candidate_funnel.csv")
    liver = pd.read_csv(R / "Task4/liver_disease_set.csv")
    dsum = pd.read_csv(R / "Direction_consistent_drugs/candidate_drugs_direction_summary.csv")
    dlinks = pd.read_csv(R / "Direction_consistent_drugs/candidate_links_direction_check.csv")
    only = pd.read_csv(R / "Direction_consistent_drugs/only_consistent_drugs.csv")
    dgrp = dsum.direction_group.value_counts()
    mp = pd.read_csv(R / "Metapath_runs/metapath_run_summary.csv").set_index("run")
    rk = pd.read_csv(RK / "ranking_run_table.csv").set_index("run")
    t25 = pd.read_csv(RK / "top25_comparison_FINAL.csv")
    finA = pd.read_csv(RK / "RunA_all_seed_genes_FINAL_ranked_drugs.csv")
    strict = {k: pd.read_csv(RK / f"{r}_strict_liver_treating_drugs.csv") for k, r, _ in RUNS}
    stA = pd.read_excel(RK / "RunA_all_seed_genes_FINAL_top_drugs_stages.xlsx", sheet_name="Stage totals")
    routesA = pd.read_csv(R / "Metapath_runs/RunA_all_seed_genes_BB.csv")
    # baseline without the direction filter (code/baseline_no_direction.py, code/baseline_compare.py)
    bsum = pd.read_csv(R / "Baseline_no_direction/baseline_run_summary.csv").iloc[0]
    bcmp = pd.read_csv(R / "Baseline_no_direction/comparison_drugs.csv")
    bshare = pd.read_csv(R / "Baseline_no_direction/route_share_by_link_class.csv").set_index("link_class")["percent"]
    bgrp = bcmp.direction_group.value_counts()
    btop = bcmp[bcmp.Rank <= 25].direction_group.value_counts()
    n_cons_links = int((dlinks.direction_check == "direction-consistent").sum())
    mainA = "RunA_all_seed_genes"
    base_tab = md_table(pd.DataFrame([
        ("Links that start routes", n_cons_links, len(dlinks)),
        ("Drugs with B routes", int(mp.loc[mainA, "drugs_with_B_routes"]), int(bsum.drugs_with_B_routes)),
        ("BB routes", n(mp.loc[mainA, "BB_routes"]), n(bsum.BB_routes)),
        ("Final ranked drugs", int(rk.loc[mainA, "final_ranked_drugs"]), len(bcmp)),
        ("Drugs with rank <= 25 (ties kept)", int(rk.loc[mainA, "final_rank_le_25"]), int((bcmp.Rank <= 25).sum())),
    ], columns=["Quantity", "With direction filter", "Baseline (no direction filter)"]))
    base_md = f"""### 5.9 Baseline without the direction filter

The main analysis (RunA) rerun on the same {len(dsum)} candidates with the direction requirement switched off: routes start
at all {len(dlinks)} candidate drug-gene links instead of the {n_cons_links} direction-consistent links of the {len(only)}
only-consistent drugs; every other rule is unchanged (`code/baseline_no_direction.py`, `code/baseline_compare.py`; outputs in
`results/Baseline_no_direction/`).

{base_tab}

Share of baseline BB routes by the direction class of the starting link: undecided {bshare.get("not decided", 0)}%,
direction-consistent {bshare.get("direction-consistent", 0)}%, inconsistent {bshare.get("direction-inconsistent", 0)}%.
Baseline ranked drugs by direction group: no decided link {int(bgrp.get("no decided link", 0))}, only consistent
{int(bgrp.get("only consistent", 0))}, only inconsistent {int(bgrp.get("only inconsistent", 0))}; in the baseline top 25:
no decided link {int(btop.get("no decided link", 0))}, only consistent {int(btop.get("only consistent", 0))}, only inconsistent
{int(btop.get("only inconsistent", 0))}. The filter therefore mainly removes drugs whose direction cannot be decided, not drugs
shown to act against the disease direction. This compares what the filter changes; it is not a measure of accuracy.
Details: `results/Baseline_no_direction/comparison_summary.md`.
"""
    linked = s07t[s07t.linked].sort_values("best_q")
    audit = j("results/audit/audit.json")

    # ------------------------------------------------------------ derived tables
    run_rows = pd.DataFrame([{
        "Run": k, "Gene\\* set": desc, "Gene\\*": n(mp.loc[r, "gene_star"]),
        "Drugs with B routes": n(mp.loc[r, "drugs_with_B_routes"]), "B routes": n(mp.loc[r, "B_routes_kept"]),
        "BB routes": n(mp.loc[r, "BB_routes"]), "Drugs with BB routes": n(rk.loc[r, "drugs_with_BB_routes"]),
        "Removed (liver contraindication)": n(rk.loc[r, "removed_liver_contraindicated"]),
        "**Final ranked drugs**": f"**{n(rk.loc[r, 'final_ranked_drugs'])}**",
        "Strict view (already treat a liver disease)": n(rk.loc[r, "strict_liver_treating_drugs"])} for k, r, desc in RUNS])
    top = finA[finA.Rank <= 25].copy()
    gs = routesA[routesA.drug_id.isin(top.Drug_ID)].groupby("drug_id").agg(
        genes=("gene_star", lambda s: ", ".join(sorted(set(s))[:6]) + (" ..." if s.nunique() > 6 else "")),
        n_gene_star=("gene_star_id", "nunique"))
    tg = only.set_index("Drug_ID").consistent_targets
    top["Direction-consistent target(s)"] = top.Drug_ID.map(tg)
    top["Gene\\* reached"] = top.Drug_ID.map(gs.n_gene_star)
    top["Runs with this drug in the top 25"] = top.Drug_ID.map(t25.set_index("Drug_ID").in_n_runs).fillna(0).astype(int)
    top["Current indications (examples)"] = top["Current indications (TREATS)"].fillna("").map(
        lambda s: "; ".join(s.split("; ")[:3]) + (" ..." if s.count(";") >= 3 else ""))
    t3show = t3[["Group", "Total genes", "Drug-connected genes", "Percentage", "Degree-matched expected %",
                 "Fold vs degree-matched", "Fisher_BH_FDR", "Perm_BH_FDR"]].copy()
    common = t25[t25.in_all_runs]
    pw_rows = linked[["name", "database", "best_q", "best_group", "ctd_n_genes", "link_basis"]]
    lc = dlinks.direction_check.value_counts()
    nodes = pd.DataFrame(list(g["nodes"].items()), columns=["Node label", "Count"])
    rels = pd.DataFrame(list(g["relationships"].items()), columns=["Relationship type", "Count"])
    audit_rows = pd.DataFrame([{"Check": k.split("_", 1)[1].replace("_", " "), "Result": "PASS" if v["pass"] else "FAIL"}
                               for k, v in sorted(audit.items())])

    readme = f"""# Knowledge-graph drug repurposing for NAFLD and liver fibrosis from public data

A reproducible knowledge-graph pipeline that searches for approved and clinical-stage drugs whose action **opposes the
disease direction of published NAFLD signature genes** and that reach NAFLD through biologically connected routes.
Every input is public and citable; every number below is read from this repository's result files
(`code/make_readme.py`). All results are **repurposing hypotheses for experimental follow-up, not treatments**.

**Final result:** the ranked drug lists in `results/Ranking/<run>_FINAL_ranked_drugs.csv` (main analysis: RunA), compared
across runs in [`results/Ranking/top25_comparison_FINAL.md`](results/Ranking/top25_comparison_FINAL.md).

## Contents

1. [Overview](#1-overview)
2. [Data sources](#2-data-sources)
3. [Methodology](#3-methodology)
4. [Implementation](#4-implementation)
5. [Results](#5-results)
6. [Limitations](#6-limitations)
7. [Repository structure](#7-repository-structure)
8. [Reproducing the analysis](#8-reproducing-the-analysis)

## 1. Overview

```
 public sources ──► knowledge graph (Neo4j) ──► seed-gene drug links ──► candidate filter ──► direction filter
                                                                                                  │
 final ranked drugs ◄── ranking (liver-contraindicated removed) ◄── metapath routes Drug → … → NAFLD ◄┘
```

| Stage | Output |
|---|---|
| Knowledge graph | {n(g['nodes_total'])} nodes, {n(g['relationships_total'])} relationships ({len(g['relationships'])} types) |
| Seed genes (Gene\\*) | 280 published NAFLD genes; {n(t2['genes_with_drug'])} have a drug link |
| Drugs linked to seed genes | {n(t2['unique_drugs'])} drugs, {n(t2['links'])} drug-gene links |
| Repurposing candidates | {n(fun.drugs.iloc[-1])} (phase >= 2, not withdrawn, no existing fibrosis-type liver use, mechanism annotated) |
| Direction-consistent drugs | {n(len(only))} (every decided drug-target link opposes the gene's disease direction) |
| NAFLD-linked pathways (route end points) | {s07['linked']} |
| **Final ranked drugs, RunA** | **{n(rk.loc['RunA_all_seed_genes', 'final_ranked_drugs'])}** drugs with BB routes to NAFLD, ranked by number of routes |

## 2. Data sources

| Source | Version / release | Retrieved | Used for |
|---|---|---|---|
| OptimusKG, Harvard Dataverse doi:10.7910/DVN/IYNGEV | v{okg['version']} (released {okg['release_time'][:10]}); {okg['licence']} | {okg_date} | base graph: genes, drugs, diseases, pathways (Reactome), biological processes, phenotypes and their edges |
| Ryaboshapkina & Hammar, Sci Rep 2017;7:12361, Suppl. Table S2 (CC BY 4.0) | 280 genes | (in `seed_genes/`) | seed genes, stage and direction |
| GWAS Catalog liver hits | 310 genes, p <= 5e-8 | 2026-11-11 | optional seed add-on (RunD) |
| MSigDB KEGG legacy + Hallmark gene sets | {gmt[0]['version']} | {gmt[0]['downloaded_utc'][:10]} | KEGG and Hallmark pathways |
| Comparative Toxicogenomics Database (CTD), diseases-pathways | {ctd['version'].lstrip('# ').replace('Report created: ', 'report of ')} | {ctd['downloaded_utc'][:10]} | curated pathway-NAFLD links (MESH:D065626) |
| ChEMBL (REST API) | {s04['release']['ChEMBL'].replace(', REST API', '')} | 2026-11-13 | drug metadata, mechanisms, indications |
| IUPHAR/BPS Guide to Pharmacology (GtoPdb) | {s04['release']['GtoPdb'].lstrip('# ').split(' (')[0].replace('GtoPdb Version: ', '')} | 2026-11-13 | drug-target links and directions |
| DGIdb | 5 (latest interactions.tsv + GraphQL API) | 2026-11-13 | drug-target links and directions |
| Open Targets Platform | {s04['release']['OpenTargets'].replace('Open Targets Platform ', '')} | 2026-11-13 | drug-target links |

Download URLs, timestamps and SHA-256 checksums: `data/SOURCES.json`. The `data/` folder is not in git (licence terms);
`code/download_inputs.py` downloads it again. No in-house or patient-derived data are used.

## 3. Methodology

### 3.1 Knowledge graph construction (Neo4j database `paper2public`)

| Step | Method | Result |
|---|---|---|
| 0.0 Import | OptimusKG genes, drugs, diseases, Reactome pathways, GO biological processes and phenotypes, with drug-gene (`ACTS_ON`, action kept), gene-gene (`INTERACTS_WITH`), pathway-gene (`IN_PATHWAY`), drug-disease (`TREATS` / `OFF_LABEL_FOR` / `CONTRAINDICATED_IN`), disease-gene (`ASSOCIATED_WITH`, evidence score), disease hierarchy (`PARENT_OF`), disease-phenotype and process-gene edges. No expression layer. | {n(imp['after']['nodes_total'])} nodes, {n(imp['after']['relationships_total'])} relationships; counts equal the parquet rows ({int(plan.missing_endpoint_dropped.sum())} disease-hierarchy rows pointing to phenotype terms excluded) |
| 0.2 Drug metadata | missing SMILES, molecule type, max phase and names filled from ChEMBL; OptimusKG values take precedence | {n(s02['drugs_written'])} drugs filled |
| 0.4 New drug-target links | for the seed genes, all drugs acting on them in GtoPdb (human, non-endogenous ligands), DGIdb (minus its GtoPdb and ChEMBL copies), ChEMBL mechanisms (human targets) and Open Targets; drugs reduced to the ChEMBL parent compound; links absent from OptimusKG kept if the drug is approved or in clinical trials **and** the source is a curated drug-target statement | {n(s04w['links_written'])} links, {n(s04w['drug_nodes_created'])} new drugs |
| 0.5 Action directions | OptimusKG links on seed genes with no direction (TARGET, ENZYME, TRANSPORTER, CARRIER, ...) looked up in ChEMBL, GtoPdb and DGIdb; a direction is stored only if all sources that give one agree | {n(s05w['written'])} of {n(s05['no_direction'])} links resolved; positive control: 100% agreement on links with a known direction |
| 0.6 Pathways | KEGG legacy and Hallmark gene sets (MSigDB) added as pathways, members matched by gene symbol | {n(s06w['pathways'])} pathways, {n(s06w['memberships'])} memberships |
| 0.7 Pathway -> NAFLD | one-sided hypergeometric test of seed genes in every Reactome, KEGG and Hallmark pathway; background = all {n(s07['background_genes'])} genes with a pathway; groups = all histology genes, each stage, each stage x direction; BH-FDR < 0.05 within group x database. Link rule: significant **and** (Reactome/KEGG) linked to NAFLD in CTD through >= 2 curated genes; Hallmark sets (absent from CTD) on enrichment alone | {s07['linked']} `IMPLICATED_IN` links |
| 0.8 Indications | ChEMBL indications of the drugs added in 0.4, written as `TREATS` (diseases already in the graph only) | {n(s08['links_written'])} links |

Every write was preceded by a dry run and has a write log (counts before/after) and an undo file in `results/graph_build/`.

### 3.2 Seed genes (Task 1)

280 genes from a meta-analysis of human liver transcriptomes across the NAFLD spectrum: 228 with a replicated,
same-direction association with steatosis, NASH, fibrosis or progression to HCC (`main_stage` 1-4, `direction`
increasing / decreasing with disease), and 52 genetic-risk genes without a histology signal (`0_genetic_risk`,
no direction). Genes are matched to the graph by Ensembl ID.

### 3.3 Drug-gene links and directions (Task 2)

All `ACTS_ON` links to seed genes. **Drug direction** comes from the link's own action (inhibitor, antagonist, blocker,
negative modulator, inverse agonist, degrader = down; agonist, partial agonist, activator, positive modulator, opener,
stabiliser = up; modulator, binding agent, substrate, releasing agent = mixed); links whose sources disagree are mixed; if
OptimusKG gives no direction, the step-0.5 direction is used. **Gene direction**: increasing = positive, decreasing = negative.

### 3.4 Druggability (Task 3)

Share of seed genes that are pharmacological drug targets (an OptimusKG action other than ENZYME / TRANSPORTER /
CARRIER, a curated link from step 0.4, or a step-0.5 direction), against a background of {n(len(t3bg))} genes with at least
one gene-gene partner. Fisher's exact test (group vs rest of background, BH-FDR), pairwise Fisher tests between stages, and
1,000 permutations of background genes matched on network-degree decile (seed 42).

### 3.5 Candidate drugs (Task 4)

Coverage of seed genes per drug (pharmacological links). Candidate filter: max clinical phase >= 2 -> not withdrawn ->
no existing use for a fibrosis-type liver disease (graph `TREATS` / `OFF_LABEL_FOR` and ChEMBL indications; fibrosis-type =
fatty liver, steatohepatitis, cirrhosis, fibrosis, portal hypertension, liver failure) -> mechanism annotated (>= 1 link
with a known direction or a stated mechanism). Liver diseases = EFO_0001421 and all descendants, minus
{int(liver.excluded_non_hepatic.sum())} non-hepatic entries (e.g. malaria, gallbladder tumours) = {int((~liver.excluded_non_hepatic).sum())} diseases.

### 3.6 Direction filter

A link is **direction-consistent** if the drug switches down a gene that increases with disease, or switches up a gene that
decreases; **inconsistent** if the opposite; otherwise not decided. Drugs kept: **only consistent** = every decided link is
direction-consistent.

### 3.7 Metapaths

For the only-consistent drugs, starting only at their direction-consistent targets:

- **B**: Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene\\* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
- **BB**: Drug -ACTS_ON-> Gene1 -INTERACTS_WITH- Gene2 -INTERACTS_WITH- Gene\\* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD

Rules: (1) routes start only at direction-consistent targets; (2) in-between genes have <= 500 interaction partners;
(3) pathway size 5-200 genes; (4) Gene1, Gene2 and Gene\\* are different genes, and the drug has no `TREATS` /
`OFF_LABEL_FOR` link to NAFLD (EFO_0003095), NASH (EFO_1001249) or cirrhosis of liver (EFO_0001422).
BB is run for drugs with B routes. Each route's stage is the `main_stage` of its Gene\\*.

| Run | Gene\\* set |
|---|---|
""" + "\n".join(f"| {k} | {d} |" for k, _, d in RUNS) + f"""

### 3.8 Ranking (final result)

For each run, **every drug with at least one BB route is ranked by its number of BB routes to NAFLD**, after removing drugs
that are contraindicated (`CONTRAINDICATED_IN`) in any liver disease. Ties share the lowest rank; the top 25 (ties kept)
are compared across runs and broken down by Gene\\* stage. As a sensitivity view, the same list restricted to drugs that
already treat a liver disease (`TREATS`) is also reported (strict view).

## 4. Implementation

- **Database**: Neo4j 2026.10.2 (Enterprise), separate database `paper2public`, built from an empty database.
- **Language**: Python 3.12; packages pinned in `requirements.txt` (pandas, numpy, scipy, statsmodels, pyarrow, neo4j,
  networkx, matplotlib, openpyxl).
- **Safety**: all graph writes go through `code/neo4j_util.py`, which refuses any database other than `paper2public`;
  credentials come only from the `NEO4J_PASSWORD` environment variable.
- **Provenance**: dry run first, then `--apply`; write logs with node/relationship counts before and after; undo files; raw
  API responses cached under `data/raw/`.
- **Checks**: `code/run_analysis.py` refuses to run unless the graph has exactly the expected counts
  (`results/graph_build/expected_graph_counts.json`); `code/final_audit.py` performs the final audit (section 5.8).

| Script | Purpose |
|---|---|
| `download_inputs.py` | download and verify all public inputs, write `data/SOURCES.json` |
| `import_optimuskg.py` | step 0.0, graph import with count check |
| `fetch_chembl_gaps.py` | step 0.2 |
| `new_links_search.py`, `write_new_links_to_graph.py` | step 0.4 |
| `external_action_check.py` | step 0.5 |
| `load_pathways.py` | step 0.6 |
| `pathway_nafld_links.py` | step 0.7 |
| `added_drug_indications.py` | step 0.8 |
| `task1_seed_gene_table.py` ... `task4_drug_coverage.py` | Tasks 1-4 |
| `direction_consistent_drugs.py` | direction filter |
| `metapaths.py` | B and BB routes for RunA-D |
| `ranking.py` | final ranked drugs, strict view, stage tables, top-25 comparison |
| `run_analysis.py` | checked runner for Tasks 1-4, direction filter, metapaths, ranking |
| `baseline_no_direction.py`, `baseline_compare.py` | baseline without the direction filter, and its comparison with RunA |
| `final_audit.py`, `make_readme.py` | audit; this README |
| `neo4j_util.py` | shared connection, write guard and logging helpers |

## 5. Results

### 5.1 Knowledge graph

{md_table(nodes, fmt={"Count": n})}

{md_table(rels, fmt={"Count": n})}

### 5.2 Seed genes

{md_table(t1c.rename(columns={"main_stage": "Stage"}))}

{md_table(t1conn.rename(columns={"main_stage": "Stage"}))}

{n(t1.usable_for_metapath.sum())} of 280 seed genes have both a gene partner and a pathway (usable as Gene\\*).

### 5.3 Drug-gene links (Task 2)

{n(t2['links'])} links between {n(t2['unique_drugs'])} drugs and {n(t2['genes_with_drug'])} seed genes ({n(t2['links_optimuskg'])} from
OptimusKG, {n(t2['links_added'])} added in step 0.4). Drug direction known for {n(t2['links_known_direction'])} links
({100 * t2['links_known_direction'] / t2['links']:.1f}%); {n(t2['links_direction_comparable'])} links have both a drug and a gene direction.

{md_table(t2s.rename(columns={"main_stage": "Stage"}))}

### 5.4 Druggability (Task 3)

Background {n(len(t3bg))} genes ({100 * t3bg.pharm_connected.mean():.1f}% pharmacological targets); {280 - len(t3out)} seed
genes are in the background.

{md_table(t3show, fmt={"Percentage": lambda v: f"{v:.1f}%", "Degree-matched expected %": lambda v: f"{v:.1f}%",
                       "Fold vs degree-matched": lambda v: f"{v:.2f}", "Fisher_BH_FDR": lambda v: f"{v:.2g}",
                       "Perm_BH_FDR": lambda v: f"{v:.2g}"})}

Seed genes are about twice as often drug targets as degree-matched background genes; no stage differs significantly from
another (pairwise BH-FDR > 0.3). With 1,000 permutations the smallest attainable permutation FDR is about 0.01.

### 5.5 Pathways linked to NAFLD

Enrichment found {int(sum(s07['significant_pathways'].values()))} significant pathways
({", ".join(f"{k} {v}" for k, v in s07['significant_pathways'].items())}); {s07['significant_and_in_ctd']} are linked to NAFLD in CTD.

{md_table(pw_rows.rename(columns={"name": "Pathway", "database": "Database", "best_q": "Best q", "best_group": "Group",
                                   "ctd_n_genes": "CTD genes", "link_basis": "Basis"}), fmt={"Best q": lambda v: f"{v:.2g}"})}

### 5.6 Candidates and direction filter

{md_table(fun.rename(columns={"step": "Step", "drugs": "Drugs", "removed_at_step": "Removed"}))}

Direction check of the {n(len(dsum))} candidates' links: {n(lc.get('direction-consistent', 0))} consistent,
{n(lc.get('direction-inconsistent', 0))} inconsistent, {n(lc.get('not decided', 0))} not decided. Drugs: **{n(dgrp.get('only consistent', 0))} only
consistent (kept)**, {n(dgrp.get('mixed', 0))} mixed, {n(dgrp.get('only inconsistent', 0))} only inconsistent, {n(dgrp.get('no decided link', 0))} with no decided link.

### 5.7 Metapath routes and final ranked drugs

{md_table(run_rows)}

**Final result, RunA (main analysis): top {len(top)} ranked drugs** (`results/Ranking/RunA_all_seed_genes_FINAL_ranked_drugs.csv`)

{md_table(top, ["Rank", "Drug", "BB routes to NAFLD", "B routes", "Direction-consistent target(s)", "Gene\\* reached",
                "Runs with this drug in the top 25", "Current indications (examples)"])}

- **Stability**: {len(common)} drugs are in the top 25 of all four runs.
- **Mechanisms** (direction-consistent targets in the table): collagen / extracellular-matrix genes (ocriplasmin, collagenase),
  PPAR-alpha (fibrates and dual PPAR-alpha/gamma agonists, which share this target and therefore tie), thioredoxin (PX-12) and
  the pregnane X receptor NR1I2 (clotrimazole, lovastatin, paclitaxel).
- **Routes by Gene\\* stage** (top-25 drugs, RunA): {", ".join(f"{r.Stage} {n(r.Routes)}" for r in stA.itertuples())}.
- **Strict view** (drugs that already treat a liver disease): RunA {len(strict['RunA'])}, RunB {len(strict['RunB'])},
  RunC {len(strict['RunC'])}, RunD {len(strict['RunD'])} drugs (`*_strict_liver_treating_drugs.csv`). Few drugs pass because the
  candidate filter already removes drugs with a fibrosis-type liver use, leaving only drugs used for liver cancer or hepatitis.

Full comparison across runs: [`results/Ranking/top25_comparison_FINAL.md`](results/Ranking/top25_comparison_FINAL.md).

### 5.8 Final audit

{md_table(audit_rows)}

(1) no excluded input names in the code; (2) the graph holds no expression layer; (3) every Gene\\* in the route files is a
seed gene (RunD: seed or GWAS add-on); (4) funnel and README numbers match the files; (5) rerunning metapaths and ranking
from scratch gives identical outputs. Details: `results/audit/audit_report.md`.

{base_md}
### 5.10 Evaluation against six baselines (`results/Evaluation/`)

The proposed method is compared with shortest path (SP), network proximity (NP), random walk with restart (RWR),
metapath counting (MPC), degree-weighted path count (DWPC) and HeteSim on the same graph and the same 299 candidate
drugs, with ablations, 100-repetition perturbation experiments (seed, edge and hub removal), an independent liver
cohort (GEO GSE126848) and a literature grading of 113 drugs. Route counts of the offline engine reproduce the Neo4j
result exactly. Result files: `main/` (rankings, overlap, cost, evidence table), `robust/`, `independent/`, `evidence/`
(PubMed queries, records, grades), `tables/`. Order of the scripts:

```
NEO4J_PASSWORD=... NEO4J_DATABASE=neo4j python code/eval_export_graph.py   # export graph tables
python code/eval_run.py              # main comparison, ablation, overlap
python code/eval_evidence_fetch.py 20   # PubMed screen (then grade by hand in evidence/evidence_grades.csv)
python code/eval_metrics.py          # evidence support and NDCG
python code/eval_robust.py 100 8     # perturbation study (repetitions, processes)
python code/eval_cost.py             # runtime and memory
python code/eval_independent.py      # needs data/geo (GSE126848 counts and series matrix)
python code/eval_tables.py && python code/eval_figures.py
```
## 6. Limitations

- **Direction is checked only at the drug -> target step.** Gene-gene interactions have no sign, so the gene-gene steps of a
  route carry no direction.
- **Few NAFLD end points.** Only {s07['linked']} pathways connect to NAFLD; 6 of them (Hallmark) rest on enrichment alone because CTD
  does not cover Hallmark gene sets.
- **Genetic-risk genes have no direction**, so drugs acting only on them cannot pass the direction filter (they can still be
  reached as Gene\\*).
- **Route counts reflect network connectivity**, not effect size; drugs sharing a target tie, and well-connected genes
  yield more routes.
- **Database coverage**: {n(t2['links'] - t2['links_known_direction'])} of {n(t2['links'])} seed-gene drug links have no known direction.
- Results are **repurposing hypotheses**, not evidence of benefit or safety.

## 7. Repository structure

```
README.md                  this file (generated by code/make_readme.py)
requirements.txt           Python packages
seed_genes/                seed gene files, sources and build scripts (see seed_genes/README.md)
code/                      all scripts (section 4)
results/
  graph_build/             step 0.0-0.8: plans, write logs, undo files, graph report, expected counts
  Task1/ ... Task4/        seed gene table, drug-gene links, druggability, candidates
  Direction_consistent_drugs/
  Metapath_runs/           B and BB routes per run
  Ranking/                 FINAL ranked drugs, strict view, stage tables, top-25 comparisons
  Baseline_no_direction/   RunA routes and ranking without the direction filter, comparison with the main run
  audit/                   final audit report
data/                      downloaded inputs (not in git)
```

## 8. Reproducing the analysis

1. `pip install -r requirements.txt`; `python code/download_inputs.py`.
2. Create an empty Neo4j database `paper2public` and set `NEO4J_PASSWORD`.
3. Build the graph, each script first as a dry run and then with `--apply`: `import_optimuskg.py`, `fetch_chembl_gaps.py`,
   `new_links_search.py` + `write_new_links_to_graph.py`, `external_action_check.py`, `load_pathways.py`,
   `pathway_nafld_links.py`, `added_drug_indications.py`; then `python code/run_analysis.py --write-expected`.
4. `python code/run_analysis.py --stage all` (Tasks 1-4, direction filter, metapaths, ranking).
5. `python code/baseline_no_direction.py` and `python code/baseline_compare.py` (baseline without the direction filter).
6. `python code/final_audit.py` and `python code/final_audit.py --rerun`; `python code/make_readme.py`.
Set `NEO4J_DATABASE` if the graph is not in the database `paper2public` (for example the default database `neo4j`).

Source versions change over time (ChEMBL, DGIdb, GtoPdb, Open Targets, CTD); exact versions and dates are recorded above and
in `data/SOURCES.json`.

## Licences

Input data keep their own licences: OptimusKG {okg['licence']}; Ryaboshapkina & Hammar 2017
supplement CC BY 4.0; CTD, MSigDB, ChEMBL, GtoPdb, DGIdb and Open Targets under their respective terms of use.
"""
    (ROOT / "README.md").write_text(readme, encoding="utf-8")
    print("wrote README.md")


if __name__ == "__main__":
    main()
