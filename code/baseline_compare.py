"""Compare the baseline without the direction requirement with the main analysis (RunA).

Reads the files written by code/baseline_no_direction.py and by the main run; no database needed.
Writes results/Baseline_no_direction/comparison_summary.md and comparison_drugs.csv.

    python code/baseline_compare.py
"""
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / "results/Baseline_no_direction"
RUN = "RunA_all_seed_genes"

def md_table(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(out)


main_fin = pd.read_csv(ROOT / f"results/Ranking/{RUN}_FINAL_ranked_drugs.csv")
main_bb = pd.read_csv(ROOT / f"results/Metapath_runs/{RUN}_BB.csv", low_memory=False)
base_fin = pd.read_csv(B / f"{RUN}_FINAL_ranked_drugs.csv")
base_bb = pd.read_csv(B / f"{RUN}_BB.csv", low_memory=False)
dsum = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_drugs_direction_summary.csv")
links = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_links_direction_check.csv")
base_sum = pd.read_csv(B / "baseline_run_summary.csv").iloc[0]
main_sum = pd.read_csv(ROOT / "results/Metapath_runs/metapath_run_summary.csv").set_index("run").loc[RUN]

group = dsum.set_index("Drug_ID").direction_group
base_fin["direction_group"] = base_fin.Drug_ID.map(group)
base_fin["retained_by_filter"] = base_fin.Drug_ID.isin(main_fin.Drug_ID)
main_rank = main_fin.set_index("Drug_ID").Rank
base_fin["main_rank"] = base_fin.Drug_ID.map(main_rank)
base_fin["main_BB_routes"] = base_fin.Drug_ID.map(main_fin.set_index("Drug_ID")["BB routes to NAFLD"])

# route share by the direction class of the link that starts the route (drug, Gene1)
lk = links[["Drug_ID", "Gene_ID", "direction_check"]].drop_duplicates()
bb = base_bb.merge(lk, left_on=["drug_id", "gene1_id"], right_on=["Drug_ID", "Gene_ID"], how="left")
route_share = (bb.direction_check.value_counts(dropna=False) / len(bb) * 100).round(1)
route_share.rename("percent").rename_axis("link_class").to_csv(B / "route_share_by_link_class.csv")

top25_base = set(base_fin[base_fin.Rank <= 25].Drug_ID)
top25_main = set(main_fin[main_fin.Rank <= 25].Drug_ID)
grp_all = base_fin.direction_group.value_counts()
grp_top = base_fin[base_fin.Rank <= 25].direction_group.value_counts()
common = base_fin[base_fin.retained_by_filter]
rho = spearmanr(common["BB routes to NAFLD"], common["main_BB_routes"])
rank_rho = spearmanr(common.Rank, common.main_rank)

# which starting targets carry the baseline routes
tg = base_bb.groupby("gene1").size().sort_values(ascending=False)
top_targets = [(g, int(n), round(100 * n / len(base_bb), 1)) for g, n in tg.head(8).items()]
# baseline top-10 drugs
top10 = base_fin.sort_values(["Rank", "Drug"]).head(10)[["Rank", "Drug", "BB routes to NAFLD", "direction_group"]]
# how many of the retained 38 sit in the baseline top 25
ret_in_top25 = int(base_fin[base_fin.retained_by_filter & (base_fin.Rank <= 25)].shape[0])
# best-ranked drug that the filter removes
removed = base_fin[~base_fin.retained_by_filter].sort_values(["Rank", "Drug"])
first_removed_rank = int(removed.Rank.min())
n_removed_top25 = int((removed.Rank <= 25).sum())
removed_top = removed.head(8)[["Rank", "Drug", "BB routes to NAFLD", "direction_group"]]
inc_top = base_fin[(base_fin.direction_group == "only inconsistent") & (base_fin.Rank <= 25)]

rows = [
    ("Candidate drugs", 321, 321),
    ("Drug-target links used to start routes", int((links.direction_check == "direction-consistent").sum()), len(links)),
    ("Drugs with B routes", int(main_sum.drugs_with_B_routes), int(base_sum.drugs_with_B_routes)),
    ("B routes", int(main_sum.B_routes_kept), int(base_sum.B_routes_kept)),
    ("BB routes", int(main_sum.BB_routes), int(base_sum.BB_routes)),
    ("Drugs removed for liver contraindication", len(main_bb.drug.unique()) - len(main_fin), len(base_bb.drug.unique()) - len(base_fin)),
    ("Final ranked drugs", len(main_fin), len(base_fin)),
    ("Drugs with rank <= 25 (ties kept)", int((main_fin.Rank <= 25).sum()), int((base_fin.Rank <= 25).sum())),
]
tab = pd.DataFrame(rows, columns=["quantity", "main (direction filter)", "baseline (no direction filter)"])
base_fin.sort_values(["Rank", "Drug"]).to_csv(B / "comparison_drugs.csv", index=False)

md = [f"# Baseline without the direction filter vs main analysis ({RUN})", "",
      "Same 321 candidates, same rules; the baseline starts routes at all 409 links, the main run only at the direction-consistent links of the 115 'only consistent' drugs.", "",
      md_table(tab), "",
      "## Ranked drugs of the baseline by direction group", "",
      "| group | ranked drugs | in top 25 (ties kept) |", "|---|---|---|"]
for g in ["only consistent", "only inconsistent", "no decided link"]:
    md.append(f"| {g} | {int(grp_all.get(g, 0))} | {int(grp_top.get(g, 0))} |")
md += ["", f"Baseline BB routes by the direction class of the starting link (%): {route_share.to_dict()}",
       f"Main-run top-25 drugs also in the baseline top 25: {len(top25_main & top25_base)} of {len(top25_main)}",
       f"Retained (direction-filter) drugs in the baseline top 25: {ret_in_top25} of {len(common)}",
       f"Drugs removed by the filter that would rank in the baseline top 25: {n_removed_top25} (best rank {first_removed_rank})",
       f"Retained drugs: Spearman rho of BB routes baseline vs main = {rho.statistic:.3f} (p={rho.pvalue:.2g}); of ranks = {rank_rho.statistic:.3f}",
       f"Top starting targets in the baseline (gene, routes, % of routes): {top_targets}", "",
       "## Baseline top 10", "", top10.to_string(index=False), "",
       "## Best-ranked drugs the filter removes", "", removed_top.to_string(index=False), ""]
(B / "comparison_summary.md").write_text("\n".join(md), encoding="utf-8")
print("\n".join(md))
