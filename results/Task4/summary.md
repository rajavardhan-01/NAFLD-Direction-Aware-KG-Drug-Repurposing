# Task 4: drug coverage and repurposing candidates

Seed genes: 280 (Task 1). Drug links: Task 2 (pharmacological-target links for the main analysis). Generated 2026-11-15T08:59:41+00:00.

## Candidate funnel

| step | drugs | removed_at_step |
|---|---|---|
| Drugs with >= 1 link to a seed gene (any action) | 1242 | 0 |
| ... with >= 1 pharmacological-target link (coverage table) | 814 | 428 |
| ... max phase >= 2 | 580 | 234 |
| ... and not withdrawn | 562 | 18 |
| ... and no existing fibrosis-type liver use = repurposing candidates | 516 | 46 |
| ... and mechanism annotated = main list | 321 | 195 |

- Liver diseases: `EFO_0001421` + PARENT_OF descendants = 231 terms, minus 70
  non-hepatic-tissue entries = 161 liver diseases; 21 are fibrosis-type (`liver_disease_set.csv`).
- Drugs hitting seed genes with any liver TREATS/OFF_LABEL_FOR use: 134;
  with a fibrosis-type liver use (excluded): 48;
  kept although they have another liver use: 85.
- Candidates with a liver contraindication (kept here, removed in the ranking step): 48.
- Drugs added in step 0.4 among the coverage table: 31; approved with no phase number (treated as phase 4):
  21; no phase number at all (cannot confirm phase >= 2, excluded): 0.
  Their existing uses come from ChEMBL indications (ChEMBL_37, 2026-05-01, retrieved 2026-11-13).

## Main list: first 20 mechanism-annotated candidates (order = histology seed genes hit, then target fraction)

| mechanism_order | Drug | max_phase | n_seed_genes | n_histology_genes | 3_fibrosis | n_down_actions | n_up_actions | target_genes |
|---|---|---|---|---|---|---|---|---|
| 1 | OCRIPLASMIN | 4.0 | 6 | 6 | 5 | 6 | 0 | COL1A2;COL3A1;COL5A1;COL6A2;LAMA2;LAMA3 |
| 2 | COLLAGENASE CLOSTRIDIUM HISTOLYTICUM | 4.0 | 4 | 4 | 4 | 4 | 0 | COL1A2;COL3A1;COL5A1;COL6A2 |
| 3 | PX-12 | 2.0 | 2 | 2 | 0 | 1 | 0 | TXN;TXNRD1 |
| 4 | FENOFIBRIC ACID | 4.0 | 3 | 2 | 0 | 0 | 1 | NR1I2;PPARA;PPARG |
| 5 | NITRIC OXIDE | 4.0 | 2 | 2 | 1 | 0 | 1 | GUCY1A1;MT1A |
| 6 | DESLANOSIDE | 4.0 | 2 | 2 | 1 | 2 | 0 | ATP1A1;FXYD2 |
| 7 | DIGITOXIN | 4.0 | 2 | 2 | 1 | 2 | 0 | ATP1A1;FXYD2 |
| 8 | ACETYLDIGITOXIN | 4.0 | 2 | 2 | 1 | 2 | 0 | ATP1A1;FXYD2 |
| 9 | LANATOSIDE C | 4.0 | 2 | 2 | 1 | 2 | 0 | ATP1A1;FXYD2 |
| 10 | ISTAROXIME | 2.0 | 2 | 2 | 1 | 2 | 0 | ATP1A1;FXYD2 |
| 11 | PACLITAXEL | 4.0 | 2 | 2 | 1 | 0 | 1 | NR1I2;TOP2A |
| 12 | TELMISARTAN | 4.0 | 2 | 1 | 0 | 1 | 0 | AGTR1;PPARG |
| 13 | NAVEGLITAZAR | 2.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 14 | MURAGLITAZAR | 3.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 15 | RAGAGLITAZAR | 2.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 16 | TESAGLITAZAR | 3.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 17 | REGLITAZAR | 2.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 18 | MK-0767 | 3.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 19 | IMIGLITAZAR | 3.0 | 2 | 1 | 0 | 0 | 2 | PPARA;PPARG |
| 20 | CANDESARTAN CILEXETIL | 4.0 | 1 | 1 | 0 | 1 | 0 | AGTR1 |

This is candidate generation, not a ranking of drug quality.

Outputs: `drug_coverage.csv` (every drug with a seed-gene link), `top50_by_coverage.csv`, `candidate_funnel.csv`,
`repurposing_candidates.csv`, `repurposing_candidates_mechanism_annotated.csv` (main list), `liver_disease_set.csv`,
`network_nodes.csv` / `network_edges.csv` / `network_top30_candidates.png`. Gene annotation: `top_enrichment_pathway`
(most significant step 0.7 pathway per seed gene). Script: `code/task4_drug_coverage.py`.
