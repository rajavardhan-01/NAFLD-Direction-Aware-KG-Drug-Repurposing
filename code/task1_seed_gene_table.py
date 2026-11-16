"""Task 1: seed gene table.

The 280 seed genes (seed_genes/seed_genes_main_Ryaboshapkina2017.csv) with their stage and
direction from the seed file, and their connections in this run's graph (paper2public, read-only):
  n_gene_partners     distinct INTERACTS_WITH partners
  n_pathways          IN_PATHWAY pathways (all), and by source (Reactome / KEGG_legacy / Hallmark)
  n_nafld_pathways    pathways with IMPLICATED_IN -> NAFLD (EFO_0003095)
  n_drugs             distinct drugs with ACTS_ON (any), of which from OptimusKG / added in step 0.4
  usable_for_metapath >= 1 gene partner and >= 1 pathway (recomputed on this graph)
No gene is dropped.

Outputs: results/Task1/seed_gene_table.csv, counts_by_stage_direction.csv, summary.md

    python code/task1_seed_gene_table.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import NAFLD, ROOT, driver, now, read_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
OUT = ROOT / "results/Task1"
STAGE_ORDER = ["0_genetic_risk", "1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]
DIR_ORDER = ["increasing", "decreasing", "mixed", "not_applicable"]

Q = """
UNWIND $ids AS id
MATCH (g:Gene {node_id: id})
CALL (g) { OPTIONAL MATCH (g)-[:INTERACTS_WITH]-(o:Gene) WHERE o <> g RETURN count(DISTINCT o) AS partners }
CALL (g) { OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)
           RETURN count(DISTINCT p) AS pathways,
                  count(DISTINCT CASE WHEN p.node_id STARTS WITH 'REACT:' THEN p END) AS reactome,
                  count(DISTINCT CASE WHEN p.source = 'KEGG_legacy' THEN p END) AS kegg,
                  count(DISTINCT CASE WHEN p.source = 'Hallmark' THEN p END) AS hallmark }
CALL (g) { OPTIONAL MATCH (g)-[:IN_PATHWAY]->(p:Pathway)-[:IMPLICATED_IN]->(:Disease {node_id: $nafld})
           RETURN count(DISTINCT p) AS nafld_pw }
CALL (g) { OPTIONAL MATCH (d:Drug)-[r:ACTS_ON]->(g)
           RETURN count(DISTINCT d) AS drugs,
                  count(DISTINCT CASE WHEN r.layer IS NULL THEN d END) AS drugs_okg,
                  count(DISTINCT CASE WHEN r.layer IS NOT NULL THEN d END) AS drugs_added }
RETURN id AS ensembl_id, g.symbol AS graph_symbol, partners, pathways, reactome, kegg, hallmark, nafld_pw,
       drugs, drugs_okg, drugs_added
"""


def md(df):
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * df.shape[1]]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.values]
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    assert len(seeds) == 280 and seeds.ensembl_id.is_unique
    drv = driver()
    with read_session(drv) as s:
        g = pd.DataFrame(s.run(Q, ids=list(seeds.ensembl_id), nafld=NAFLD).data())
    drv.close()
    missing = set(seeds.ensembl_id) - set(g.ensembl_id)
    if missing:
        raise SystemExit(f"{len(missing)} seed genes not in graph: {sorted(missing)[:5]}")
    g = g.rename(columns={"partners": "n_gene_partners", "pathways": "n_pathways", "reactome": "n_pathways_reactome",
                          "kegg": "n_pathways_kegg", "hallmark": "n_pathways_hallmark", "nafld_pw": "n_nafld_pathways",
                          "drugs": "n_drugs", "drugs_okg": "n_drugs_optimuskg", "drugs_added": "n_drugs_added"})
    t = seeds[["gene", "ensembl_id", "entrez", "main_stage", "all_stages", "evidence", "direction", "source",
               "also_in_GWAS_Catalog"]].merge(g, on="ensembl_id", how="left")
    t["usable_for_metapath"] = (t.n_gene_partners > 0) & (t.n_pathways > 0)
    t["symbol_matches_graph"] = t.gene.eq(t.graph_symbol)
    assert len(t) == 280
    t.to_csv(OUT / "seed_gene_table.csv", index=False)

    ct = pd.crosstab(t.main_stage, t.direction, margins=True, margins_name="all")
    ct = ct.reindex(index=STAGE_ORDER + ["all"], columns=[c for c in DIR_ORDER if c in ct.columns] + ["all"], fill_value=0)
    ct.to_csv(OUT / "counts_by_stage_direction.csv")
    conn = t.groupby("main_stage").agg(genes=("gene", "size"), usable_for_metapath=("usable_for_metapath", "sum"),
                                       with_partner=("n_gene_partners", lambda s: int((s > 0).sum())),
                                       with_pathway=("n_pathways", lambda s: int((s > 0).sum())),
                                       with_nafld_pathway=("n_nafld_pathways", lambda s: int((s > 0).sum())),
                                       with_drug=("n_drugs", lambda s: int((s > 0).sum()))).reindex(STAGE_ORDER)
    conn.loc["all"] = conn.sum()
    conn.to_csv(OUT / "connections_by_stage.csv")
    all_stage = t.all_stages.str.split(";").explode().str.strip().value_counts().reindex(STAGE_ORDER, fill_value=0)

    summary = f"""# Task 1: seed gene table

{len(t)} seed genes from Ryaboshapkina & Hammar, Sci Rep 2017;7:12361 (Supplementary Table S2), all found in the graph
(symbol differs from the graph symbol for {int((~t.symbol_matches_graph).sum())} genes; matched by Ensembl ID).
Generated {now()} from database paper2public.

## Genes by main stage and direction

{md(ct.reset_index())}

Genes listing each stage in `all_stages` (a gene can list several): {", ".join(f"{k} {v}" for k, v in all_stage.items())}.

## Connections in the graph

{md(conn.reset_index().astype({c: int for c in conn.columns}))}

- usable_for_metapath = at least one INTERACTS_WITH partner and at least one pathway.
- Drug links: {int((t.n_drugs > 0).sum())} genes have at least one drug ({int((t.n_drugs_optimuskg > 0).sum())} via OptimusKG,
  {int((t.n_drugs_added > 0).sum())} via links added in step 0.4).
- NAFLD-linked pathways (step 0.7): {int((t.n_nafld_pathways > 0).sum())} seed genes are members of at least one.

Output: `seed_gene_table.csv` (one row per seed gene), `counts_by_stage_direction.csv`, `connections_by_stage.csv`.
Script: `code/task1_seed_gene_table.py`.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
