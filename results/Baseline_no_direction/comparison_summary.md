# Baseline without the direction filter vs main analysis (RunA_all_seed_genes)

Same 321 candidates, same rules; the baseline starts routes at all 409 links, the main run only at the direction-consistent links of the 115 'only consistent' drugs.

| quantity | main (direction filter) | baseline (no direction filter) |
|---|---|---|
| Candidate drugs | 321 | 321 |
| Drug-target links used to start routes | 128 | 409 |
| Drugs with B routes | 44 | 130 |
| B routes | 166 | 487 |
| BB routes | 5025 | 21324 |
| Drugs removed for liver contraindication | 6 | 14 |
| Final ranked drugs | 38 | 116 |
| Drugs with rank <= 25 (ties kept) | 27 | 38 |

## Ranked drugs of the baseline by direction group

| group | ranked drugs | in top 25 (ties kept) |
|---|---|---|
| only consistent | 40 | 13 |
| only inconsistent | 2 | 0 |
| no decided link | 74 | 25 |

Baseline BB routes by the direction class of the starting link (%): {'not decided': 74.4, 'direction-consistent': 24.8, 'direction-inconsistent': 0.8}
Main-run top-25 drugs also in the baseline top 25: 13 of 27
Retained (direction-filter) drugs in the baseline top 25: 13 of 38
Drugs removed by the filter that would rank in the baseline top 25: 25 (best rank 3)
Retained drugs: Spearman rho of BB routes baseline vs main = 0.799 (p=1.9e-09); of ranks = 0.799
Top starting targets in the baseline (gene, routes, % of routes): [('PPARG', 8760, 41.1), ('TNF', 3332, 15.6), ('PPARA', 2400, 11.3), ('TGFB1', 1015, 4.8), ('NR1I2', 846, 4.0), ('COL1A2', 402, 1.9), ('ATP1A1', 375, 1.8), ('TOP2A', 372, 1.7)]

## Baseline top 10

 Rank                                 Drug  BB routes to NAFLD direction_group
    1                          OCRIPLASMIN                 633 only consistent
    2 COLLAGENASE CLOSTRIDIUM HISTOLYTICUM                 595 only consistent
    3                           MESALAMINE                 401 no decided link
    4                         INDOMETHACIN                 352 no decided link
    5                          ALEGLITAZAR                 315 only consistent
    5                             GW501516                 315 no decided link
    5                          IMIGLITAZAR                 315 only consistent
    5                         INDEGLITAZAR                 315 only consistent
    5                              MK-0767                 315 only consistent
    5                         MURAGLITAZAR                 315 only consistent

## Best-ranked drugs the filter removes

 Rank          Drug  BB routes to NAFLD direction_group
    3    MESALAMINE                 401 no decided link
    4  INDOMETHACIN                 352 no decided link
    5      GW501516                 315 no decided link
   17   CANNABIDIOL                 224 no decided link
   18  ARHALOFENATE                 219 no decided link
   18     ATX08-001                 219 no decided link
   18 BALAGLITAZONE                 219 no decided link
   18   BALSALAZIDE                 219 no decided link
