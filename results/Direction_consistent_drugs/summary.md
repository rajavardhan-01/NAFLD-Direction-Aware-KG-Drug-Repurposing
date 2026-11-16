# Direction filter

Input: 321 mechanism-annotated repurposing candidates (Task 4) and their 409 drug-gene links (Task 2).

## Links
- direction-consistent: 128
- direction-inconsistent: 43
- not decided: 238 (96 with mixed/unknown drug direction; 142 with a known drug
  direction on a 0_genetic_risk gene, which has no direction)

## Drugs (321)

| group | drugs |
|---|---|
| only consistent | 115 |
| mixed | 0 |
| only inconsistent | 43 |
| no decided link | 163 |
| total | 321 |

**Kept for the metapaths: 115 "only consistent" drugs** (every decided link is direction-consistent).
12 of them have a liver contraindication (removed later in the ranking step).

Seed genes most often hit by a consistent link (only-consistent drugs): TOP2A 27, PPARA 23, CA12 14, TYMS 12, NR1I2 7, FXYD2 5, ATP1A1 5, TXNRD1 4, ITGB5 3, COL3A1 2.

Rule: consistent = drug pushes DOWN a gene that increases with disease, or UP a gene that decreases. Directions are checked only
at the drug -> target step. This is a mechanism-consistency filter, not evidence that a drug helps.

Outputs: `only_consistent_drugs.csv/.xlsx` (metapath input), `direction_consistent_drugs.csv` (only consistent + mixed),
`direction_consistent_links.csv`, `candidate_drugs_direction_summary.csv`, `candidate_links_direction_check.csv`.
Script: `code/direction_consistent_drugs.py`.
