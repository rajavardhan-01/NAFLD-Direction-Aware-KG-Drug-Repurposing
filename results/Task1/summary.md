# Task 1: seed gene table

280 seed genes from Ryaboshapkina & Hammar, Sci Rep 2017;7:12361 (Supplementary Table S2), all found in the graph
(symbol differs from the graph symbol for 0 genes; matched by Ensembl ID).
Generated 2026-11-16T11:32:16+00:00 from database paper2public.

## Genes by main stage and direction

| main_stage | increasing | decreasing | not_applicable | all |
|---|---|---|---|---|
| 0_genetic_risk | 0 | 0 | 52 | 52 |
| 1_steatosis | 10 | 3 | 0 | 13 |
| 2_NASH | 56 | 27 | 0 | 83 |
| 3_fibrosis | 93 | 29 | 0 | 122 |
| 4_HCC | 1 | 9 | 0 | 10 |
| all | 160 | 68 | 52 | 280 |

Genes listing each stage in `all_stages` (a gene can list several): 0_genetic_risk 52, 1_steatosis 13, 2_NASH 89, 3_fibrosis 146, 4_HCC 41.

## Connections in the graph

| main_stage | genes | usable_for_metapath | with_partner | with_pathway | with_nafld_pathway | with_drug |
|---|---|---|---|---|---|---|
| 0_genetic_risk | 52 | 51 | 52 | 51 | 14 | 32 |
| 1_steatosis | 13 | 12 | 13 | 12 | 2 | 4 |
| 2_NASH | 83 | 68 | 81 | 69 | 15 | 28 |
| 3_fibrosis | 122 | 93 | 120 | 93 | 30 | 35 |
| 4_HCC | 10 | 9 | 10 | 9 | 4 | 6 |
| all | 280 | 233 | 276 | 234 | 65 | 105 |

- usable_for_metapath = at least one INTERACTS_WITH partner and at least one pathway.
- Drug links: 105 genes have at least one drug (103 via OptimusKG,
  29 via links added in step 0.4).
- NAFLD-linked pathways (step 0.7): 65 seed genes are members of at least one.

Output: `seed_gene_table.csv` (one row per seed gene), `counts_by_stage_direction.csv`, `connections_by_stage.csv`.
Script: `code/task1_seed_gene_table.py`.
