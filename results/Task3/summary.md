# Task 3: druggability of the seed genes

Background: **18,453 genes** with at least one INTERACTS_WITH partner in this graph
(pharmacological target 17.5%, any drug 18.8%).
276 of 280 seed genes are in the background; 4 have no partner and are left out
(`seed_genes_outside_background.csv`). Drug flags for seed genes match Task 2. Generated 2026-11-15T08:42:18+00:00.

## A. Pharmacological target (main definition)

| Group | Total genes | Drug-connected genes | Percentage | Degree-matched expected % | Fold vs degree-matched | Odds ratio | Fisher_BH_FDR | Perm_BH_FDR |
|---|---|---|---|---|---|---|---|---|
| All seed genes | 276 | 100 | 36.2 | 18.8 | 1.92 | 2.73 | 1.23e-12 | 0.00999 |
| Histology genes (stages 1-4) | 224 | 69 | 30.8 | 18.7 | 1.65 | 2.12 | 6.36e-06 | 0.00999 |
| 0_genetic_risk | 52 | 31 | 59.6 | 19.7 | 3.03 | 7.01 | 9.42e-11 | 0.00999 |
| 1_steatosis | 13 | 4 | 30.8 | 20.1 | 1.53 | 2.09 | 0.328 | 0.611 |
| 2_NASH | 81 | 26 | 32.1 | 18.8 | 1.7 | 2.24 | 0.00521 | 0.00999 |
| 3_fibrosis | 120 | 33 | 27.5 | 17.6 | 1.56 | 1.79 | 0.014 | 0.0214 |
| 4_HCC | 10 | 6 | 60 | 22.3 | 2.69 | 7.07 | 0.00679 | 0.0524 |
| 1_steatosis_increasing | 10 | 3 | 30 | 17.7 | 1.7 | 2.02 | 0.455 | 0.611 |
| 1_steatosis_decreasing | 3 | 1 | 33.3 | 26.9 | 1.24 | 2.35 | 0.47 | 1 |
| 2_NASH_increasing | 55 | 15 | 27.3 | 19.5 | 1.4 | 1.77 | 0.11 | 0.27 |
| 2_NASH_decreasing | 26 | 11 | 42.3 | 17.6 | 2.41 | 3.46 | 0.00679 | 0.00999 |
| 3_fibrosis_increasing | 91 | 24 | 26.4 | 18 | 1.47 | 1.69 | 0.0611 | 0.111 |
| 3_fibrosis_decreasing | 29 | 9 | 31 | 16.4 | 1.9 | 2.12 | 0.111 | 0.111 |
| 4_HCC_increasing | 1 | 0 | 0 | 15.6 | 0 | 0 | 1 | 1 |
| 4_HCC_decreasing | 9 | 6 | 66.7 | 22.8 | 2.92 | 9.43 | 0.00521 | 0.00999 |

Pairwise stage comparisons (Fisher, BH-FDR): see `pairwise_stage_comparisons_pharmacological.csv`.

| Comparison | %A | %B | Odds ratio | Fisher_p | BH_FDR |
|---|---|---|---|---|---|
| 1_steatosis vs 2_NASH | 30.8 | 32.1 | 0.941 | 1 | 1 |
| 1_steatosis vs 3_fibrosis | 30.8 | 27.5 | 1.17 | 0.755 | 0.905 |
| 1_steatosis vs 4_HCC | 30.8 | 60 | 0.314 | 0.222 | 0.443 |
| 2_NASH vs 3_fibrosis | 32.1 | 27.5 | 1.25 | 0.529 | 0.794 |
| 2_NASH vs 4_HCC | 32.1 | 60 | 0.32 | 0.157 | 0.443 |
| 3_fibrosis vs 4_HCC | 27.5 | 60 | 0.256 | 0.0649 | 0.389 |

## B. Any drug link (secondary)

| Group | Total genes | Drug-connected genes | Percentage | Degree-matched expected % | Fold vs degree-matched | Odds ratio | Fisher_BH_FDR | Perm_BH_FDR |
|---|---|---|---|---|---|---|---|---|
| All seed genes | 276 | 105 | 38 | 20.2 | 1.89 | 2.7 | 1.04e-12 | 0.00749 |
| Histology genes (stages 1-4) | 224 | 73 | 32.6 | 20 | 1.63 | 2.11 | 3.58e-06 | 0.00749 |
| 0_genetic_risk | 52 | 32 | 61.5 | 20.9 | 2.94 | 6.95 | 9.72e-11 | 0.00749 |
| 1_steatosis | 13 | 4 | 30.8 | 21.3 | 1.45 | 1.92 | 0.355 | 0.701 |
| 2_NASH | 81 | 28 | 34.6 | 20.1 | 1.72 | 2.29 | 0.00317 | 0.00749 |
| 3_fibrosis | 120 | 35 | 29.2 | 18.9 | 1.54 | 1.78 | 0.0124 | 0.0214 |
| 4_HCC | 10 | 6 | 60 | 23.4 | 2.57 | 6.47 | 0.00995 | 0.0566 |
| 1_steatosis_increasing | 10 | 3 | 30 | 19.1 | 1.57 | 1.85 | 0.475 | 0.701 |
| 1_steatosis_decreasing | 3 | 1 | 33.3 | 27.8 | 1.2 | 2.15 | 0.498 | 1 |
| 2_NASH_increasing | 55 | 16 | 29.1 | 20.6 | 1.41 | 1.77 | 0.0864 | 0.237 |
| 2_NASH_decreasing | 26 | 12 | 46.1 | 19.1 | 2.42 | 3.7 | 0.00405 | 0.012 |
| 3_fibrosis_increasing | 91 | 24 | 26.4 | 19.3 | 1.37 | 1.55 | 0.108 | 0.195 |
| 3_fibrosis_decreasing | 29 | 11 | 37.9 | 17.8 | 2.13 | 2.64 | 0.025 | 0.0375 |
| 4_HCC_increasing | 1 | 0 | 0 | 17.1 | 0 | 0 | 1 | 1 |
| 4_HCC_decreasing | 9 | 6 | 66.7 | 24 | 2.78 | 8.63 | 0.00552 | 0.02 |

Definitions: pharmacological target = an OptimusKG action other than ENZYME / TRANSPORTER / CARRIER, a curated link added in step 0.4,
or an OptimusKG link carrying a step 0.5 direction. Fisher's exact test compares a group with the rest of the background.
The permutation test draws 1,000 random background gene sets matched on INTERACTS_WITH degree decile (seed 42); well-studied genes
have more partners and more drugs, so the degree-matched comparison is the fairer one.

Outputs: `druggability_summary_pharmacological.csv`, `druggability_summary_any.csv`, `pairwise_stage_comparisons_*.csv`,
`background_genes_druggability.csv`, `background_degree_bins.csv`, `druggability_barplot.png`. Script: `code/task3_druggability.py`.
