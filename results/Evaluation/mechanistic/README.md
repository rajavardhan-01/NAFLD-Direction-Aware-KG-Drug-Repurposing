# Independent pharmacological and mechanistic evaluation of the prioritised drugs (steps 1-6)

Faculty request: export the 116 conventional-MPC drugs, identify the 38 retained and verify the 78 excluded; extract the
drug-target pairs and actions; collect independent experimental publications; grade drug-target-disease evidence with the
M0-M4 framework; compare retained / excluded / MPC / annotation-only groups; report confidence intervals, contradictory
evidence and target-family sensitivity.

**Important:** the M0-M4 levels below are *our working definitions* (the faculty's own wording was not available). Replace
them in `code/mech_step4_grade.py` if the faculty defines the levels differently, then rerun steps 4-6.

## Run order (all read-only on the graph)

```
NEO4J_PASSWORD=... NEO4J_DATABASE=neo4j python code/mech_step1_2_export.py   # steps 1-2
python code/mech_step3_collect.py        # step 3: PubMed search, about 25 minutes, needs internet
python code/mech_step4_grade.py          # step 4: M0-M4 grades from the records
python code/mech_step5_6_compare.py      # steps 5-6: groups, intervals, contradictions, family analysis
```
Step 5-6 reads `step4_manual_review.csv` (the hand review of the automatic contradiction flags) and rewrites
`step4_pair_grades.csv` with the confirmed flags, so always run step 4 before step 5-6.

## Steps 1-2: groups and pairs
- `step1_drug_groups.csv`: the 116 MPC drugs (offline list equals the Neo4j list), 38 retained, 78 excluded with reasons:
  74 have no decided direction on any seed link; 2 (olmesartan, lorlatinib) are only inconsistent; 2 (daunorubicin,
  celecoxib) are consistent, but their consistent targets (TOP2A, CA12) have no B route and BB is only searched after a
  B route. The annotation-only group (decidability-only filter, drugs with routes) has 39 drugs, 38 of them retained.
- `step2_drug_target_pairs.csv`: 532 pairs (all ACTS_ON targets) with action label, direction, source; 163 pairs are on
  seed genes (116 drugs, 40 genes). Steps 3-6 use these 163 seed-gene pairs, the pairs that start the routes.

## Step 3: independent publications
`step3_pubmed_records.json`: for each of the 163 pairs two PubMed queries (drug terms AND target terms; the same plus a
liver-disease term), up to 40 records each. "Independent" = primary research articles from PubMed, not the pharmacology
databases that built the graph and not the source signature paper. Retrieved 2026-11-22.

## Step 4: M0-M4 framework (rule based, record-level flags in `step4_record_flags.csv`)
| Grade | Meaning |
|---|---|
| M0 | no primary publication mentions both the drug and the target |
| M1 | mention only, no experiment on the drug at the target |
| M2 | experimental engagement: binding, enzyme/receptor/reporter assay, IC50/EC50/Ki |
| M3 | experimental pharmacological action in cells/animals/humans, in a sentence naming drug, target and an action that matches the annotated direction (any action if the direction is unknown or mixed) |
| M4 | M3 plus the same publication names NAFLD / NASH / steatosis / liver fibrosis |

Drug-level grade = best M over the drug's seed-gene pairs. Reviews, editorials, comments, letters, meta-analyses and
retracted papers are not counted.

Validation by hand: a first, looser version gave 93 of 163 pairs M3 and was rejected after an audit (genetics papers,
clinical-trial mentions). The sentence-level version was audited again on 14 random M3 pairs: 11 supported, 3 not
(simlukafusp alfa-FAP is weak, paclitaxel-TOP2A cites another molecule, clofibrate-CYP2E1 is indirect), so precision of
M3 is roughly 80%. The 42 automatic contradiction flags were all read by hand (`step4_manual_review.csv`): 6 confirmed
(eritoran-TLR4 in two drug entries, bardoxolone methyl-PPARG, indomethacin-PPARA, ocriplasmin-LAMA2, ocriplasmin-LAMA3),
36 rejected.

## Steps 5-6: results (drug unit, 95% Wilson intervals)
| Group | n | M>=2 | M>=3 | M4 | contradictory (confirmed) |
|---|---|---|---|---|---|
| retained (direction filter) | 38 | 26 (68%) | 25 (66%, 50-79) | 10 (26%) | 1 (3%) |
| excluded | 78 | 57 (73%) | 48 (62%, 50-72) | 12 (15%) | 4 (5%) |
| conventional MPC | 116 | 83 (72%) | 73 (63%, 54-71) | 22 (19%) | 5 (4%) |
| annotation-only | 39 | 26 (67%) | 25 (64%, 48-77) | 10 (26%) | 1 (3%) |

Retained minus excluded (bootstrap 95% interval, Fisher p): M>=3 +0.04 (-0.15 to +0.22, p=0.69); M4 +0.11 (-0.05 to
+0.28, p=0.21). Annotation-only vs retained: no difference (M>=3 -0.02, -0.23 to +0.19, p=1.0; M4 -0.01, -0.19 to +0.20, p=1.0), see `step6_difference_tests.csv`.

Target-family sensitivity (`step6_family_sensitivity.csv`): with one value per target family, M>=3 is 12 of 14 families
(86%) for retained and 19 of 22 (86%) for excluded. Without the PPARG, PPARA and TNF families, retained drugs fall to
M>=3 10 of 17 (59%) and M4 1 of 17 (6%); excluded drugs to 20 of 37 (54%) and 6 of 37 (16%).

**Reading:** the drugs excluded by the direction filter have mechanistic evidence as well supported as the retained drugs;
the filter does not select drugs with better-supported mechanisms, and the retained group's higher M4 share rests on the
PPAR drugs and is not statistically different. All grades come from keyword rules on PubMed abstracts (limits: abstracts
only, English-language keyword matching, single rule set, no blinded second grader), so they are screening grades, not
expert curation.

Files: `step5_drug_level.csv`, `step5_group_comparison.csv`, `step6_difference_tests.csv`, `step6_family_sensitivity.csv`,
`step6_contradictions.csv`, `step4_pair_grades.csv`, `step4_record_flags.csv`.
