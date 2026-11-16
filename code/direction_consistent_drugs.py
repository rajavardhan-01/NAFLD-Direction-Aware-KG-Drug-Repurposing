"""Direction filter: repurposing candidates whose drug action opposes the seed gene's disease direction.

Input: the Task 4 main list (mechanism-annotated repurposing candidates) and all their Task 2 drug-gene links.
A link is
  * direction-consistent   : the drug switches the gene product DOWN and the gene is 'increasing' with disease,
                             or switches it UP and the gene is 'decreasing';
  * direction-inconsistent : the opposite;
  * not decided            : the drug's direction is mixed / unknown, or the gene has no direction
                             (0_genetic_risk genes are 'not_applicable').
Drug groups: only consistent | mixed | only inconsistent | no decided link.
Kept for the metapaths: "only consistent" drugs (every decided link is direction-consistent).
This is a mechanism-consistency filter, not evidence that a drug helps. No Neo4j needed.

Outputs: results/Direction_consistent_drugs/.

    python code/direction_consistent_drugs.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/Direction_consistent_drugs"
EDGES = ROOT / "results/Task2/drug_gene_edges.csv"
MAIN = ROOT / "results/Task4/repurposing_candidates_mechanism_annotated.csv"
GROUPS = ["only consistent", "mixed", "only inconsistent", "no decided link"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    edges = pd.read_csv(EDGES, low_memory=False)
    main_ = pd.read_csv(MAIN, low_memory=False)
    e = edges[edges.Drug_ID.isin(main_.Drug_ID)].copy()
    known = e.Action_direction.isin(["down", "up"]) & e.gene_direction.isin(["positive", "negative"])
    e["direction_check"] = "not decided"
    e.loc[known & (((e.Action_direction == "down") & (e.gene_direction == "positive"))
                   | ((e.Action_direction == "up") & (e.gene_direction == "negative"))), "direction_check"] = "direction-consistent"
    e.loc[known & (((e.Action_direction == "down") & (e.gene_direction == "negative"))
                   | ((e.Action_direction == "up") & (e.gene_direction == "positive"))), "direction_check"] = "direction-inconsistent"
    e.to_csv(OUT / "candidate_links_direction_check.csv", index=False)

    g = e.groupby("Drug_ID")
    drugs = pd.DataFrame({
        "n_consistent": g.direction_check.apply(lambda s: int((s == "direction-consistent").sum())),
        "n_inconsistent": g.direction_check.apply(lambda s: int((s == "direction-inconsistent").sum())),
        "n_not_decided": g.direction_check.apply(lambda s: int((s == "not decided").sum())),
        "consistent_targets": e[e.direction_check == "direction-consistent"].groupby("Drug_ID").Gene.agg(lambda s: ";".join(sorted(set(s)))),
        "inconsistent_targets": e[e.direction_check == "direction-inconsistent"].groupby("Drug_ID").Gene.agg(lambda s: ";".join(sorted(set(s)))),
    }).reset_index()
    drugs["has_known_direction"] = (drugs.n_consistent + drugs.n_inconsistent) > 0
    drugs["direction_group"] = "no decided link"
    drugs.loc[(drugs.n_consistent > 0) & (drugs.n_inconsistent == 0), "direction_group"] = "only consistent"
    drugs.loc[(drugs.n_consistent > 0) & (drugs.n_inconsistent > 0), "direction_group"] = "mixed"
    drugs.loc[(drugs.n_consistent == 0) & (drugs.n_inconsistent > 0), "direction_group"] = "only inconsistent"
    cols = ["Drug_ID", "Drug", "max_phase", "is_approved", "added_drug_node", "liver_contraindicated",
            "has_liver_indication", "liver_indications", "n_seed_genes", "n_histology_genes"]
    drugs = drugs.merge(main_[cols], on="Drug_ID", how="left")
    assert len(drugs) == len(main_), "every main-list drug must appear once"
    drugs.to_csv(OUT / "candidate_drugs_direction_summary.csv", index=False)

    keep = drugs[drugs.n_consistent > 0].sort_values(["n_consistent", "n_histology_genes"], ascending=False)
    keep.to_csv(OUT / "direction_consistent_drugs.csv", index=False)
    e[e.Drug_ID.isin(keep.Drug_ID) & (e.direction_check == "direction-consistent")].to_csv(
        OUT / "direction_consistent_links.csv", index=False)
    only = keep[keep.direction_group == "only consistent"]
    only.to_csv(OUT / "only_consistent_drugs.csv", index=False)
    only.to_excel(OUT / "only_consistent_drugs.xlsx", index=False)

    counts = drugs.direction_group.value_counts().reindex(GROUPS, fill_value=0)
    lc = e.direction_check.value_counts()
    nd = e[e.direction_check == "not decided"]
    why = pd.Series({"drug direction mixed/unknown": int((~nd.Action_direction.isin(["down", "up"])).sum()),
                     "drug direction known, gene has no direction (0_genetic_risk)":
                         int((nd.Action_direction.isin(["down", "up"]) & nd.gene_direction.isna()).sum())})
    tg = only.consistent_targets.str.split(";").explode().value_counts()
    summary = f"""# Direction filter

Input: {len(main_)} mechanism-annotated repurposing candidates (Task 4) and their {len(e):,} drug-gene links (Task 2).

## Links
- direction-consistent: {lc.get('direction-consistent', 0):,}
- direction-inconsistent: {lc.get('direction-inconsistent', 0):,}
- not decided: {lc.get('not decided', 0):,} ({why.iloc[0]:,} with mixed/unknown drug direction; {why.iloc[1]:,} with a known drug
  direction on a 0_genetic_risk gene, which has no direction)

## Drugs ({len(drugs)})

| group | drugs |
|---|---|
""" + "\n".join(f"| {k} | {v} |" for k, v in counts.items()) + f"""
| total | {int(counts.sum())} |

**Kept for the metapaths: {len(only)} "only consistent" drugs** (every decided link is direction-consistent).
{int(only.liver_contraindicated.sum())} of them have a liver contraindication (removed later in the ranking step).

Seed genes most often hit by a consistent link (only-consistent drugs): {", ".join(f"{k} {v}" for k, v in tg.head(10).items())}.

Rule: consistent = drug pushes DOWN a gene that increases with disease, or UP a gene that decreases. Directions are checked only
at the drug -> target step. This is a mechanism-consistency filter, not evidence that a drug helps.

Outputs: `only_consistent_drugs.csv/.xlsx` (metapath input), `direction_consistent_drugs.csv` (only consistent + mixed),
`direction_consistent_links.csv`, `candidate_drugs_direction_summary.csv`, `candidate_links_direction_check.csv`.
Script: `code/direction_consistent_drugs.py`.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
