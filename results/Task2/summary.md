# Task 2: drug-gene links on the seed genes

Gene set: the 280 seed genes (Task 1). Source: database paper2public, read-only
(`MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids`), queried 2026-11-14T05:50:23+00:00.

- **1,659 links** between **1,242 drugs** and **105 of 280 seed genes**
  (1,549 OptimusKG, 110 added in step 0.4). No duplicate drug-gene pairs: True.
- **Known drug direction: 400 links (24.1%)**: down 236, up 164;
  mixed 21, unknown 1,238.
- Direction source: none 1,238, OptimusKG action 292, added link (push) 67, GtoPdb/DGIdb/ChEMBL on OptimusKG link (step 0.5) 62.
- Links with a known drug direction AND a gene direction (comparable for the direction filter): 198.
  `0_genetic_risk` genes have no direction (`not_applicable`), so their links are never comparable.
- Graph totals: 18,671 Drug nodes, 36,721 ACTS_ON links; match the step 0.4 write log: True.

## By stage of the seed gene

| main_stage | links | genes | drugs | known_direction | comparable |
|---|---|---|---|---|---|
| 0_genetic_risk | 588 | 32 | 527 | 202 | 0 |
| 1_steatosis | 24 | 4 | 24 | 18 | 18 |
| 2_NASH | 439 | 28 | 412 | 22 | 22 |
| 3_fibrosis | 302 | 35 | 285 | 91 | 91 |
| 4_HCC | 306 | 6 | 284 | 67 | 67 |

## By link origin

| Link_origin | links | genes | drugs | known_direction |
|---|---|---|---|---|
| ChEMBL/DGIdb/GtoPdb/OpenTargets (step 0.4) | 110 | 29 | 108 | 64 |
| OptimusKG | 1549 | 103 | 1151 | 336 |

## Most common actions

| Action_raw | links |
|---|---|
| TARGET | 540 |
| ENZYME | 459 |
| TRANSPORTER | 204 |
| INHIBITOR | 113 |
| AGONIST | 101 |
| ANTAGONIST | 43 |
| unknown | 42 |
| CARRIER | 30 |
| Inhibitor / Inhibition | 24 |
| ACTIVATOR | 13 |
| BINDING_AGENT | 12 |
| Agonist / Agonist | 11 |
| Antagonist / Antagonist | 11 |
| HYDROLYTIC_ENZYME | 10 |
| ANTISENSE_INHIBITOR | 7 |

Direction rule: drug direction from the link's own action; OptimusKG actions INHIBITOR, ANTAGONIST, BLOCKER, NEGATIVE_(ALLOSTERIC_)MODULATOR,
ALLOSTERIC_ANTAGONIST, INVERSE_AGONIST, DEGRADER = down; AGONIST, PARTIAL_AGONIST, ACTIVATOR, POSITIVE_(ALLOSTERIC_)MODULATOR, OPENER,
STABILISER = up; MODULATOR, BINDING_AGENT, SUBSTRATE, RELEASING_AGENT = mixed; others unknown, then filled from step 0.5 where available.
Gene direction: increasing = positive, decreasing = negative.

Output: `drug_gene_edges.csv` (1,659 rows), `summary_numbers.csv`, `links_by_stage.csv`, `links_by_origin.csv`, `graph_check.json`.
Script: `code/task2_drug_gene_edges.py`.
