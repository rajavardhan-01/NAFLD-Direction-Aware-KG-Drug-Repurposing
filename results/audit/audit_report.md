# Final audit

## 1_excluded_names_in_code: PASS (2026-11-16T11:31:57+00:00)
Searched 23 files for 32 patterns; hits: 0

## 2_graph_layers: PASS (2026-11-16T11:31:57+00:00)
Database neo4j. TRACKS_* types: none; expression-like labels: none; Gene nodes with expression-like properties: 0.
Labels: BiologicalProcess 25,754, Disease 36,753, Drug 18,671, Gene 61,306, Pathway 3,041, Phenotype 19,341
Relationship types: ACTS_ON 36,721, ASSOCIATED_WITH 9,734,774, CONTRAINDICATED_IN 11,379, HAS_PHENOTYPE 157,144, IMPLICATED_IN 12, INTERACTS_WITH 327,924, IN_PATHWAY 67,044, OFF_LABEL_FOR 827, PARENT_OF 47,750, PARTICIPATES_IN 158,410, TREATS 58,509

## 3_gene_star_in_seed_file: PASS (2026-11-16T11:31:57+00:00)
- RunA_all_seed_genes: 55 Gene* in routes, outside the allowed set: 0
- RunB_histology_genes: 42 Gene* in routes, outside the allowed set: 0
- RunC_fibrosis_genes: 25 Gene* in routes, outside the allowed set: 0
- RunD_seed_plus_GWAS: 85 Gene* in routes, outside the allowed set: 0

## 4_numbers: PASS (2026-11-15T14:08:53+00:00)
- funnel_steps_add_up: True
- funnel_monotonic: True
- direction_groups_match_main_list: True
- ranking_counts_match_files: True
- readme_equals_fresh_build_from_files: True

## 5_rerun_identical: PASS (2026-11-16T11:30:47+00:00)
18 files compared; identical: True; differing: none
