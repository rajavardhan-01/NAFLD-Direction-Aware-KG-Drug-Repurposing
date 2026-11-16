# Knowledge-graph drug repurposing for NAFLD and liver fibrosis from public data

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
| Knowledge graph | 164,866 nodes, 10,600,494 relationships (11 types) |
| Seed genes (Gene\*) | 280 published NAFLD genes; 105 have a drug link |
| Drugs linked to seed genes | 1,242 drugs, 1,659 drug-gene links |
| Repurposing candidates | 321 (phase >= 2, not withdrawn, no existing fibrosis-type liver use, mechanism annotated) |
| Direction-consistent drugs | 115 (every decided drug-target link opposes the gene's disease direction) |
| NAFLD-linked pathways (route end points) | 12 |
| **Final ranked drugs, RunA** | **38** drugs with BB routes to NAFLD, ranked by number of routes |

## 2. Data sources

| Source | Version / release | Retrieved | Used for |
|---|---|---|---|
| OptimusKG, Harvard Dataverse doi:10.7910/DVN/IYNGEV | v3.2 (released 2026-09-09); CC BY-NC-SA 4.0 | 2026-11-13 | base graph: genes, drugs, diseases, pathways (Reactome), biological processes, phenotypes and their edges |
| Ryaboshapkina & Hammar, Sci Rep 2017;7:12361, Suppl. Table S2 (CC BY 4.0) | 280 genes | (in `seed_genes/`) | seed genes, stage and direction |
| GWAS Catalog liver hits | 310 genes, p <= 5e-8 | 2026-11-11 | optional seed add-on (RunD) |
| MSigDB KEGG legacy + Hallmark gene sets | 2024.1.Hs | 2026-11-13 | KEGG and Hallmark pathways |
| Comparative Toxicogenomics Database (CTD), diseases-pathways | report of Thu Oct 29 12:36:38 EDT 2026 | 2026-11-13 | curated pathway-NAFLD links (MESH:D065626) |
| ChEMBL (REST API) | ChEMBL_37 (2026-05-01) | 2026-11-13 | drug metadata, mechanisms, indications |
| IUPHAR/BPS Guide to Pharmacology (GtoPdb) | 2026.3.1 - published: 2026-10-12 | 2026-11-13 | drug-target links and directions |
| DGIdb | 5 (latest interactions.tsv + GraphQL API) | 2026-11-13 | drug-target links and directions |
| Open Targets Platform | 26.09 | 2026-11-13 | drug-target links |

Download URLs, timestamps and SHA-256 checksums: `data/SOURCES.json`. The `data/` folder is not in git (licence terms);
`code/download_inputs.py` downloads it again. No in-house or patient-derived data are used.

## 3. Methodology

### 3.1 Knowledge graph construction (Neo4j database `paper2public`)

| Step | Method | Result |
|---|---|---|
| 0.0 Import | OptimusKG genes, drugs, diseases, Reactome pathways, GO biological processes and phenotypes, with drug-gene (`ACTS_ON`, action kept), gene-gene (`INTERACTS_WITH`), pathway-gene (`IN_PATHWAY`), drug-disease (`TREATS` / `OFF_LABEL_FOR` / `CONTRAINDICATED_IN`), disease-gene (`ASSOCIATED_WITH`, evidence score), disease hierarchy (`PARENT_OF`), disease-phenotype and process-gene edges. No expression layer. | 164,599 nodes, 10,579,970 relationships; counts equal the parquet rows (142 disease-hierarchy rows pointing to phenotype terms excluded) |
| 0.2 Drug metadata | missing SMILES, molecule type, max phase and names filled from ChEMBL; OptimusKG values take precedence | 2,696 drugs filled |
| 0.4 New drug-target links | for the seed genes, all drugs acting on them in GtoPdb (human, non-endogenous ligands), DGIdb (minus its GtoPdb and ChEMBL copies), ChEMBL mechanisms (human targets) and Open Targets; drugs reduced to the ChEMBL parent compound; links absent from OptimusKG kept if the drug is approved or in clinical trials **and** the source is a curated drug-target statement | 110 links, 31 new drugs |
| 0.5 Action directions | OptimusKG links on seed genes with no direction (TARGET, ENZYME, TRANSPORTER, CARRIER, ...) looked up in ChEMBL, GtoPdb and DGIdb; a direction is stored only if all sources that give one agree | 62 of 1,257 links resolved; positive control: 100% agreement on links with a known direction |
| 0.6 Pathways | KEGG legacy and Hallmark gene sets (MSigDB) added as pathways, members matched by gene symbol | 236 pathways, 20,067 memberships |
| 0.7 Pathway -> NAFLD | one-sided hypergeometric test of seed genes in every Reactome, KEGG and Hallmark pathway; background = all 12,420 genes with a pathway; groups = all histology genes, each stage, each stage x direction; BH-FDR < 0.05 within group x database. Link rule: significant **and** (Reactome/KEGG) linked to NAFLD in CTD through >= 2 curated genes; Hallmark sets (absent from CTD) on enrichment alone | 12 `IMPLICATED_IN` links |
| 0.8 Indications | ChEMBL indications of the drugs added in 0.4, written as `TREATS` (diseases already in the graph only) | 335 links |

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
CARRIER, a curated link from step 0.4, or a step-0.5 direction), against a background of 18,453 genes with at least
one gene-gene partner. Fisher's exact test (group vs rest of background, BH-FDR), pairwise Fisher tests between stages, and
1,000 permutations of background genes matched on network-degree decile (seed 42).

### 3.5 Candidate drugs (Task 4)

Coverage of seed genes per drug (pharmacological links). Candidate filter: max clinical phase >= 2 -> not withdrawn ->
no existing use for a fibrosis-type liver disease (graph `TREATS` / `OFF_LABEL_FOR` and ChEMBL indications; fibrosis-type =
fatty liver, steatohepatitis, cirrhosis, fibrosis, portal hypertension, liver failure) -> mechanism annotated (>= 1 link
with a known direction or a stated mechanism). Liver diseases = EFO_0001421 and all descendants, minus
70 non-hepatic entries (e.g. malaria, gallbladder tumours) = 161 diseases.

### 3.6 Direction filter

A link is **direction-consistent** if the drug switches down a gene that increases with disease, or switches up a gene that
decreases; **inconsistent** if the opposite; otherwise not decided. Drugs kept: **only consistent** = every decided link is
direction-consistent.

### 3.7 Metapaths

For the only-consistent drugs, starting only at their direction-consistent targets:

- **B**: Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene\* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
- **BB**: Drug -ACTS_ON-> Gene1 -INTERACTS_WITH- Gene2 -INTERACTS_WITH- Gene\* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD

Rules: (1) routes start only at direction-consistent targets; (2) in-between genes have <= 500 interaction partners;
(3) pathway size 5-200 genes; (4) Gene1, Gene2 and Gene\* are different genes, and the drug has no `TREATS` /
`OFF_LABEL_FOR` link to NAFLD (EFO_0003095), NASH (EFO_1001249) or cirrhosis of liver (EFO_0001422).
BB is run for drugs with B routes. Each route's stage is the `main_stage` of its Gene\*.

| Run | Gene\* set |
|---|---|
| RunA | all 280 seed genes (main analysis) |
| RunB | 228 histology genes (`main_stage` != `0_genetic_risk`) |
| RunC | all genes with `3_fibrosis` in `all_stages` |
| RunD | 280 seed genes + GWAS Catalog add-on genes |

### 3.8 Ranking (final result)

For each run, **every drug with at least one BB route is ranked by its number of BB routes to NAFLD**, after removing drugs
that are contraindicated (`CONTRAINDICATED_IN`) in any liver disease. Ties share the lowest rank; the top 25 (ties kept)
are compared across runs and broken down by Gene\* stage. As a sensitivity view, the same list restricted to drugs that
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
| `final_audit.py`, `make_readme.py` | audit; this README |
| `neo4j_util.py` | shared connection, write guard and logging helpers |

## 5. Results

### 5.1 Knowledge graph

| Node label | Count |
|---|---|
| BiologicalProcess | 25,754 |
| Disease | 36,753 |
| Drug | 18,671 |
| Gene | 61,306 |
| Pathway | 3,041 |
| Phenotype | 19,341 |

| Relationship type | Count |
|---|---|
| ACTS_ON | 36,721 |
| ASSOCIATED_WITH | 9,734,774 |
| CONTRAINDICATED_IN | 11,379 |
| HAS_PHENOTYPE | 157,144 |
| IMPLICATED_IN | 12 |
| INTERACTS_WITH | 327,924 |
| IN_PATHWAY | 67,044 |
| OFF_LABEL_FOR | 827 |
| PARENT_OF | 47,750 |
| PARTICIPATES_IN | 158,410 |
| TREATS | 58,509 |

### 5.2 Seed genes

| Stage | increasing | decreasing | not_applicable | all |
|---|---|---|---|---|
| 0_genetic_risk | 0 | 0 | 52 | 52 |
| 1_steatosis | 10 | 3 | 0 | 13 |
| 2_NASH | 56 | 27 | 0 | 83 |
| 3_fibrosis | 93 | 29 | 0 | 122 |
| 4_HCC | 1 | 9 | 0 | 10 |
| all | 160 | 68 | 52 | 280 |

| Stage | genes | usable_for_metapath | with_partner | with_pathway | with_nafld_pathway | with_drug |
|---|---|---|---|---|---|---|
| 0_genetic_risk | 52 | 51 | 52 | 51 | 14 | 32 |
| 1_steatosis | 13 | 12 | 13 | 12 | 2 | 4 |
| 2_NASH | 83 | 68 | 81 | 69 | 15 | 28 |
| 3_fibrosis | 122 | 93 | 120 | 93 | 30 | 35 |
| 4_HCC | 10 | 9 | 10 | 9 | 4 | 6 |
| all | 280 | 233 | 276 | 234 | 65 | 105 |

233 of 280 seed genes have both a gene partner and a pathway (usable as Gene\*).

### 5.3 Drug-gene links (Task 2)

1,659 links between 1,242 drugs and 105 seed genes (1,549 from
OptimusKG, 110 added in step 0.4). Drug direction known for 400 links
(24.1%); 198 links have both a drug and a gene direction.

| Stage | links | genes | drugs | known_direction | comparable |
|---|---|---|---|---|---|
| 0_genetic_risk | 588 | 32 | 527 | 202 | 0 |
| 1_steatosis | 24 | 4 | 24 | 18 | 18 |
| 2_NASH | 439 | 28 | 412 | 22 | 22 |
| 3_fibrosis | 302 | 35 | 285 | 91 | 91 |
| 4_HCC | 306 | 6 | 284 | 67 | 67 |

### 5.4 Druggability (Task 3)

Background 18,453 genes (17.5% pharmacological targets); 276 seed
genes are in the background.

| Group | Total genes | Drug-connected genes | Percentage | Degree-matched expected % | Fold vs degree-matched | Fisher_BH_FDR | Perm_BH_FDR |
|---|---|---|---|---|---|---|---|
| All seed genes | 276 | 100 | 36.2% | 18.8% | 1.92 | 1.2e-12 | 0.01 |
| Histology genes (stages 1-4) | 224 | 69 | 30.8% | 18.7% | 1.65 | 6.4e-06 | 0.01 |
| 0_genetic_risk | 52 | 31 | 59.6% | 19.7% | 3.03 | 9.4e-11 | 0.01 |
| 1_steatosis | 13 | 4 | 30.8% | 20.1% | 1.53 | 0.33 | 0.61 |
| 2_NASH | 81 | 26 | 32.1% | 18.8% | 1.70 | 0.0052 | 0.01 |
| 3_fibrosis | 120 | 33 | 27.5% | 17.6% | 1.56 | 0.014 | 0.021 |
| 4_HCC | 10 | 6 | 60.0% | 22.3% | 2.69 | 0.0068 | 0.052 |
| 1_steatosis_increasing | 10 | 3 | 30.0% | 17.7% | 1.70 | 0.46 | 0.61 |
| 1_steatosis_decreasing | 3 | 1 | 33.3% | 26.9% | 1.24 | 0.47 | 1 |
| 2_NASH_increasing | 55 | 15 | 27.3% | 19.5% | 1.40 | 0.11 | 0.27 |
| 2_NASH_decreasing | 26 | 11 | 42.3% | 17.6% | 2.41 | 0.0068 | 0.01 |
| 3_fibrosis_increasing | 91 | 24 | 26.4% | 18.0% | 1.47 | 0.061 | 0.11 |
| 3_fibrosis_decreasing | 29 | 9 | 31.0% | 16.4% | 1.90 | 0.11 | 0.11 |
| 4_HCC_increasing | 1 | 0 | 0.0% | 15.6% | 0.00 | 1 | 1 |
| 4_HCC_decreasing | 9 | 6 | 66.7% | 22.8% | 2.92 | 0.0052 | 0.01 |

Seed genes are about twice as often drug targets as degree-matched background genes; no stage differs significantly from
another (pairwise BH-FDR > 0.3). With 1,000 permutations the smallest attainable permutation FDR is about 0.01.

### 5.5 Pathways linked to NAFLD

Enrichment found 23 significant pathways
(Hallmark 6, KEGG 3, Reactome 14); 9 are linked to NAFLD in CTD.

| Pathway | Database | Best q | Group | CTD genes | Basis |
|---|---|---|---|---|---|
| Epithelial mesenchymal transition | Hallmark | 2e-06 | all_histology | 0 | enrichment only (Hallmark, no CTD coverage) |
| Focal adhesion | KEGG | 0.0013 | 3_fibrosis_increasing | 4 | enrichment + CTD |
| ECM proteoglycans | Reactome | 0.0018 | all_histology | 2 | enrichment + CTD |
| Non-integrin membrane-ECM interactions | Reactome | 0.0023 | 3_fibrosis_increasing | 3 | enrichment + CTD |
| Angiogenesis | Hallmark | 0.0043 | all_histology | 0 | enrichment only (Hallmark, no CTD coverage) |
| Syndecan interactions | Reactome | 0.0044 | 3_fibrosis_increasing | 2 | enrichment + CTD |
| Xenobiotic metabolism | Hallmark | 0.011 | 4_HCC_decreasing | 0 | enrichment only (Hallmark, no CTD coverage) |
| Extracellular matrix organization | Reactome | 0.014 | 3_fibrosis_increasing | 4 | enrichment + CTD |
| E2F targets | Hallmark | 0.021 | 3_fibrosis_increasing | 0 | enrichment only (Hallmark, no CTD coverage) |
| MTORC1 signaling | Hallmark | 0.025 | 2_NASH_increasing | 0 | enrichment only (Hallmark, no CTD coverage) |
| Chemokine signaling pathway | KEGG | 0.032 | 3_fibrosis | 4 | enrichment + CTD |
| Hypoxia | Hallmark | 0.047 | all_histology | 0 | enrichment only (Hallmark, no CTD coverage) |

### 5.6 Candidates and direction filter

| Step | Drugs | Removed |
|---|---|---|
| Drugs with >= 1 link to a seed gene (any action) | 1242 | 0 |
| ... with >= 1 pharmacological-target link (coverage table) | 814 | 428 |
| ... max phase >= 2 | 580 | 234 |
| ... and not withdrawn | 562 | 18 |
| ... and no existing fibrosis-type liver use = repurposing candidates | 516 | 46 |
| ... and mechanism annotated = main list | 321 | 195 |

Direction check of the 321 candidates' links: 128 consistent,
43 inconsistent, 238 not decided. Drugs: **115 only
consistent (kept)**, 0 mixed, 43 only inconsistent, 163 with no decided link.

### 5.7 Metapath routes and final ranked drugs

| Run | Gene\* set | Gene\* | Drugs with B routes | B routes | BB routes | Drugs with BB routes | Removed (liver contraindication) | **Final ranked drugs** | Strict view (already treat a liver disease) |
|---|---|---|---|---|---|---|---|---|---|
| RunA | all 280 seed genes (main analysis) | 280 | 44 | 166 | 5,025 | 44 | 6 | **38** | 1 |
| RunB | 228 histology genes (`main_stage` != `0_genetic_risk`) | 228 | 37 | 106 | 2,675 | 37 | 2 | **35** | 0 |
| RunC | all genes with `3_fibrosis` in `all_stages` | 146 | 27 | 83 | 1,853 | 27 | 2 | **25** | 0 |
| RunD | 280 seed genes + GWAS Catalog add-on genes | 573 | 73 | 278 | 15,290 | 73 | 9 | **64** | 10 |

**Final result, RunA (main analysis): top 27 ranked drugs** (`results/Ranking/RunA_all_seed_genes_FINAL_ranked_drugs.csv`)

| Rank | Drug | BB routes to NAFLD | B routes | Direction-consistent target(s) | Gene\* reached | Runs with this drug in the top 25 | Current indications (examples) |
|---|---|---|---|---|---|---|---|
| 1 | OCRIPLASMIN | 633 | 29 | COL1A2;COL3A1;COL5A1;COL6A2;LAMA2;LAMA3 | 34 | 4 | deep vein thrombosis; eye disease; macular degeneration ... |
| 2 | COLLAGENASE CLOSTRIDIUM HISTOLYTICUM | 595 | 29 | COL1A2;COL3A1;COL5A1;COL6A2 | 32 | 4 | Dupuytren Contracture; Palmar Fibromatosis; Penile Fibromatosis ... |
| 3 | PX-12 | 130 | 1 | TXN | 32 | 3 | pancreatic neoplasm |
| 4 | ALEGLITAZAR | 96 | 3 | PPARA | 24 | 4 | type 2 diabetes mellitus |
| 4 | AVE8134 | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | CHOLINE FENOFIBRATE | 96 | 3 | PPARA | 24 | 4 | Combined hyperlipidemia; Disorder of lipid metabolism; cardiovascular disease ... |
| 4 | CIPROFIBRATE | 96 | 3 | PPARA | 24 | 4 | cardiovascular disease |
| 4 | CLOFIBRATE | 96 | 3 | PPARA | 24 | 4 | Myocardial Ischemia; cardiovascular disease; hyperlipidemia ... |
| 4 | CP-778875 | 96 | 3 | PPARA | 24 | 4 | Disorder of lipid metabolism |
| 4 | DRF-10945 | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | GW590735 | 96 | 3 | PPARA | 24 | 4 | Disorder of lipid metabolism |
| 4 | IMIGLITAZAR | 96 | 3 | PPARA | 24 | 4 | diabetes mellitus |
| 4 | INDEGLITAZAR | 96 | 3 | PPARA | 24 | 4 | diabetes mellitus |
| 4 | LY-518674 | 96 | 3 | PPARA | 24 | 4 | hypertension; metabolic syndrome |
| 4 | MK-0767 | 96 | 3 | PPARA | 24 | 4 | diabetes mellitus; metabolic syndrome; type 2 diabetes mellitus |
| 4 | MURAGLITAZAR | 96 | 3 | PPARA | 24 | 4 | Disorder of lipid metabolism; diabetes mellitus; type 2 diabetes mellitus |
| 4 | NAVEGLITAZAR | 96 | 3 | PPARA | 24 | 4 | type 2 diabetes mellitus |
| 4 | PEMAFIBRATE | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | PIRINIXIC ACID | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | RAGAGLITAZAR | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | REGLITAZAR | 96 | 3 | PPARA | 24 | 4 |  |
| 4 | SODELGLITAZAR | 96 | 3 | PPARA | 24 | 4 | metabolic syndrome; type 2 diabetes mellitus |
| 4 | TESAGLITAZAR | 96 | 3 | PPARA | 24 | 4 | type 2 diabetes mellitus |
| 4 | ZYH7 | 96 | 3 | PPARA | 24 | 4 | Disorder of lipid metabolism; Hypertriglyceridemia |
| 25 | CLOTRIMAZOLE | 94 | 1 | NR1I2 | 25 | 1 | HIV infection; bacterial vaginosis; candidiasis ... |
| 25 | LOVASTATIN | 94 | 1 | NR1I2 | 25 | 1 | Cerebral arteriovenous malformation; Disorder of lipid metabolism; HIV infection ... |
| 25 | PACLITAXEL | 94 | 1 | NR1I2 | 25 | 1 | Central Nervous System Neoplasm; Coronary Restenosis; Fallopian Tube Carcinoma ... |

- **Stability**: 23 drugs are in the top 25 of all four runs.
- **Mechanisms** (direction-consistent targets in the table): collagen / extracellular-matrix genes (ocriplasmin, collagenase),
  PPAR-alpha (fibrates and dual PPAR-alpha/gamma agonists, which share this target and therefore tie), thioredoxin (PX-12) and
  the pregnane X receptor NR1I2 (clotrimazole, lovastatin, paclitaxel).
- **Routes by Gene\* stage** (top-25 drugs, RunA): 0_genetic_risk 1,298, 1_steatosis 16, 2_NASH 456, 3_fibrosis 1,883, 4_HCC 3.
- **Strict view** (drugs that already treat a liver disease): RunA 1, RunB 0,
  RunC 0, RunD 10 drugs (`*_strict_liver_treating_drugs.csv`). Few drugs pass because the
  candidate filter already removes drugs with a fibrosis-type liver use, leaving only drugs used for liver cancer or hepatitis.

Full comparison across runs: [`results/Ranking/top25_comparison_FINAL.md`](results/Ranking/top25_comparison_FINAL.md).

### 5.8 Final audit

| Check | Result |
|---|---|
| excluded names in code | PASS |
| graph layers | PASS |
| gene star in seed file | PASS |
| numbers | PASS |
| rerun identical | PASS |

(1) no excluded input names in the code; (2) the graph holds no expression layer; (3) every Gene\* in the route files is a
seed gene (RunD: seed or GWAS add-on); (4) funnel and README numbers match the files; (5) rerunning metapaths and ranking
from scratch gives identical outputs. Details: `results/audit/audit_report.md`.

## 6. Limitations

- **Direction is checked only at the drug -> target step.** Gene-gene interactions have no sign, so the gene-gene steps of a
  route carry no direction.
- **Few NAFLD end points.** Only 12 pathways connect to NAFLD; 6 of them (Hallmark) rest on enrichment alone because CTD
  does not cover Hallmark gene sets.
- **Genetic-risk genes have no direction**, so drugs acting only on them cannot pass the direction filter (they can still be
  reached as Gene\*).
- **Route counts reflect network connectivity**, not effect size; drugs sharing a target tie, and well-connected genes
  yield more routes.
- **Database coverage**: 1,259 of 1,659 seed-gene drug links have no known direction.
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
5. `python code/final_audit.py` and `python code/final_audit.py --rerun`; `python code/make_readme.py`.

Source versions change over time (ChEMBL, DGIdb, GtoPdb, Open Targets, CTD); exact versions and dates are recorded above and
in `data/SOURCES.json`.

## Licences

Input data keep their own licences: OptimusKG CC BY-NC-SA 4.0; Ryaboshapkina & Hammar 2017
supplement CC BY 4.0; CTD, MSigDB, ChEMBL, GtoPdb, DGIdb and Open Targets under their respective terms of use.
