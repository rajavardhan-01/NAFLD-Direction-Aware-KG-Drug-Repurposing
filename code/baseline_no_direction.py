"""Baseline without the direction requirement.

Runs the same metapath routes (B, then BB) and the same ranking as the main analysis (RunA, all 280 seed genes) on the
same 321 mechanism-annotated candidate drugs, but starts the routes at ALL of each candidate's drug-gene links to seed
genes (409 links: direction-consistent, inconsistent and not decided) instead of only the direction-consistent links of
the 115 'only consistent' drugs. Every other rule is unchanged: hub limit (<= 500 partners), pathway size (5-200 genes),
Gene1 / Gene2 / Gene* different, exclusion of drugs that already treat NAFLD, NASH or cirrhosis, and the removal of
drugs contraindicated in any liver disease in the ranking.

Outputs go to results/Baseline_no_direction/ (the main outputs are not touched). Read-only on Neo4j.

    python code/baseline_no_direction.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import metapaths as mp  # noqa: E402
import ranking as rk  # noqa: E402
from neo4j_util import ROOT, driver, now, read_session  # noqa: E402

OUT = ROOT / "results/Baseline_no_direction"
LINKS = ROOT / "results/Direction_consistent_drugs/candidate_links_direction_check.csv"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mp.OUT = OUT          # write route files here instead of results/Metapath_runs
    rk.OUT = OUT          # write ranking files here instead of results/Ranking
    rk.RUNS_DIR = OUT
    links = pd.read_csv(LINKS)
    drugs = links[["Drug", "Drug_ID"]].drop_duplicates()
    pairs = (links[["Drug", "Drug_ID", "Gene_ID"]].drop_duplicates()
             .rename(columns={"Drug": "drug", "Drug_ID": "drug_id", "Gene_ID": "gene_id"}).to_dict("records"))
    name, genes = mp.gene_sets()["RunA"]
    print(f"baseline without direction filter: {len(drugs)} drugs, {len(pairs)} drug-target pairs, {len(genes)} Gene*")
    liver = pd.read_csv(rk.LIVER_SET)
    liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
    drv = driver()
    with read_session(drv) as ses:
        known = {r["id"]: r["n"] for r in ses.run(
            """UNWIND $ids AS id MATCH (d:Drug {node_id: id})
               OPTIONAL MATCH (d)-[:TREATS|OFF_LABEL_FOR]->(x:Disease) WHERE x.node_id IN $liver
               RETURN id, count(x) AS n""", ids=list(drugs.Drug_ID), liver=mp.KNOWN_LIVER)}
        res = mp.run_one(ses, name, genes, drugs, pairs, known)
        info, F, T = rk.rank_run(ses, name, liver_ids)
    drv.close()
    pd.DataFrame([{"generated": now(), **res}]).to_csv(OUT / "baseline_run_summary.csv", index=False)
    print(info)


if __name__ == "__main__":
    main()
