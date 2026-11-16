"""Task 2: every drug-gene (ACTS_ON) link on the 280 seed genes, with action, source and direction.

Source: database paper2public only (read-only). Two kinds of ACTS_ON links:
  * OptimusKG links (action in `action`, per-source assertions, sources_direct, ...)
  * links added in step 0.4 (layer = 'added curated drug-target links'; action labels in `action_labels`,
    direction in `push`, databases in `sources`)
Direction logic:
  * Drug direction = the link's own action (OptimusKG action mapped to down/up/mixed/unknown, demoted to
    'mixed' when the source assertions point both ways; for added links the stored `push`).
    If an OptimusKG link gives none, the `external_action_direction` stored on the same link (step 0.5) is used.
  * Gene direction (seed file `direction`): increasing -> positive, decreasing -> negative;
    mixed / not_applicable -> no gene direction.

Outputs: results/Task2/drug_gene_edges.csv, summary_numbers.csv, graph_check.json, summary.md

    python code/task2_drug_gene_edges.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, driver, now, read_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
OUT = ROOT / "results/Task2"
ADDED_LAYER = "added curated drug-target links"
NEW_LINKS_LOG = ROOT / "results/graph_build/step0_4_new_links/graph_write/write_log.json"

DOWN = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE_MODULATOR", "NEGATIVE_ALLOSTERIC_MODULATOR",
        "ALLOSTERIC_ANTAGONIST", "INVERSE_AGONIST", "DEGRADER"}
UP = {"AGONIST", "PARTIAL_AGONIST", "ACTIVATOR", "POSITIVE_MODULATOR", "POSITIVE_ALLOSTERIC_MODULATOR",
      "OPENER", "STABILISER"}
MIXED = {"MODULATOR", "BINDING_AGENT", "SUBSTRATE", "RELEASING_AGENT"}
PUSH = {"down": "down", "up": "up", "conflicting": "mixed", "unknown": "unknown"}
GENE_DIR = {"increasing": "positive", "decreasing": "negative"}


def action_direction(a):
    a = str(a).upper()
    return "down" if a in DOWN else "up" if a in UP else "mixed" if a in MIXED else "unknown"


def edge_direction(action, assertions):
    d = action_direction(action)
    if isinstance(assertions, str) and assertions:
        if {"down", "up"} <= {action_direction(x.split(":", 1)[-1]) for x in assertions.split(";")}:
            return "mixed"
    return d


def blank(x):
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == ""


def joined(v):
    return ";".join(map(str, v)) if isinstance(v, list) else v


QUERY = """
MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids
RETURN d.node_id AS Drug_ID, coalesce(d.recovered_name, d.name) AS Drug, d.molecule_type AS molecule_type,
       d.max_phase AS max_phase, d.is_approved AS is_approved, d.has_been_withdrawn AS withdrawn,
       d.source AS drug_node_source, d.drug_stage AS drug_stage,
       g.ensembl_id AS Gene_ID, g.symbol AS graph_symbol, type(r) AS Relationship, r.layer AS layer,
       r.action AS okg_action, r.relation_assertions AS relation_assertions, r.relation_conflict AS relation_conflict,
       r.sources_direct AS sources_direct, r.sources_indirect AS sources_indirect,
       r.evidence_level AS evidence_level, r.mechanisms_of_action AS okg_mechanism, r.source_urls AS source_urls,
       r.action_labels AS ext_action_labels, r.push AS ext_push, r.sources AS ext_sources,
       r.source_detail AS ext_source_detail, r.evidence AS ext_evidence, r.mechanism AS ext_mechanism,
       r.release AS ext_release, r.added AS ext_added, r.n_sources AS n_sources,
       r.external_action_direction AS external_action_direction, r.external_action_label AS external_action_label,
       r.external_action_sources AS external_action_sources
"""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    drv = driver()
    with read_session(drv) as s:
        edges = pd.DataFrame(s.run(QUERY, ids=list(seeds.ensembl_id)).data())
        totals = s.run("""MATCH (d:Drug) WITH count(d) AS drugs MATCH ()-[r:ACTS_ON]->(:Gene)
            RETURN drugs, count(r) AS acts_on, sum(CASE WHEN r.layer = $l THEN 1 ELSE 0 END) AS added""",
                       l=ADDED_LAYER).single()
        degree = pd.DataFrame(s.run("""MATCH (d:Drug)-[r:ACTS_ON]->(:Gene) WHERE d.node_id IN $ids
            RETURN d.node_id AS Drug_ID, count(r) AS drug_total_targets""",
                                    ids=list(edges.Drug_ID.unique())).data())
    drv.close()
    for c in ["sources_direct", "sources_indirect", "source_urls", "okg_mechanism"]:
        edges[c] = edges[c].map(joined)

    added = edges.layer.eq(ADDED_LAYER)
    edges["Link_origin"] = added.map({True: "ChEMBL/DGIdb/GtoPdb/OpenTargets (step 0.4)", False: "OptimusKG"})
    edges["Action_raw"] = [(lab if not blank(lab) else "unknown") if a else (act if not blank(act) else "unknown")
                           for a, lab, act in zip(added, edges.ext_action_labels, edges.okg_action)]
    edges["Action_direction"] = [PUSH.get(str(p), "unknown") if a else edge_direction(act, ass)
                                 for a, p, act, ass in zip(added, edges.ext_push, edges.okg_action, edges.relation_assertions)]
    fill = (~added) & edges.Action_direction.eq("unknown") & edges.external_action_direction.isin(["down", "up"])
    edges.loc[fill, "Action_direction"] = edges.loc[fill, "external_action_direction"]
    edges["Direction_source"] = "none"
    edges.loc[added & edges.Action_direction.ne("unknown"), "Direction_source"] = "added link (push)"
    edges.loc[(~added) & edges.Action_direction.ne("unknown") & ~fill, "Direction_source"] = "OptimusKG action"
    edges.loc[fill, "Direction_source"] = "GtoPdb/DGIdb/ChEMBL on OptimusKG link (step 0.5)"
    edges["Evidence_source"] = ["External/" + str(src).replace(";", "+") if a else
                                "OptimusKG/" + ("" if blank(sd) else str(sd).replace(";", "+"))
                                for a, src, sd in zip(added, edges.ext_sources, edges.sources_direct)]
    edges["mechanism_of_action"] = edges.okg_mechanism.where(~added, edges.ext_mechanism)
    edges["evidence_level"] = edges.evidence_level.where(~added, edges.ext_evidence)
    edges["Database_ID"] = edges.Drug_ID + "|" + edges.Gene_ID

    ann = seeds[["ensembl_id", "gene", "main_stage", "all_stages", "direction"]].rename(
        columns={"ensembl_id": "Gene_ID", "gene": "Gene", "direction": "seed_direction"})
    df = edges.merge(ann, on="Gene_ID", how="left").merge(degree, on="Drug_ID", how="left")
    assert len(df) == len(edges)
    df["gene_direction"] = df.seed_direction.map(GENE_DIR)
    df["direction_comparable"] = df.Action_direction.isin(["up", "down"]) & df.gene_direction.notna()
    cols = ["Drug", "Drug_ID", "Gene", "Gene_ID", "Relationship", "Action_raw", "Action_direction", "Direction_source",
            "Evidence_source", "Database_ID", "Link_origin", "evidence_level", "n_sources", "relation_assertions",
            "relation_conflict", "mechanism_of_action", "sources_indirect", "source_urls", "ext_source_detail",
            "ext_release", "ext_added", "external_action_direction", "external_action_label", "external_action_sources",
            "molecule_type", "max_phase", "is_approved", "withdrawn", "drug_stage", "drug_total_targets",
            "main_stage", "all_stages", "seed_direction", "gene_direction", "direction_comparable"]
    df = df[cols].sort_values(["Gene", "Drug"])
    df.to_csv(OUT / "drug_gene_edges.csv", index=False)

    log = json.loads(NEW_LINKS_LOG.read_text(encoding="utf-8"))
    gcheck = {"source_used": "Neo4j paper2public (read-only)", "queried": now(),
              "graph_drug_nodes": totals["drugs"], "graph_acts_on_total": totals["acts_on"],
              "graph_acts_on_added": totals["added"],
              "write_log_acts_on_after_step0_4": log["after"]["relationships"]["ACTS_ON"],
              "write_log_drugs_after_step0_4": log["after"]["nodes"]["Drug"],
              "edges_on_seed_genes": len(df), "edges_optimuskg": int((df.Link_origin == "OptimusKG").sum()),
              "edges_added": int((df.Link_origin != "OptimusKG").sum()),
              "duplicate_drug_gene_pairs": int(df.duplicated(["Drug_ID", "Gene_ID"]).sum()),
              "edges_with_external_action_props": int(df.external_action_direction.notna().sum())}
    gcheck["graph_matches_write_log"] = (gcheck["graph_acts_on_total"] == gcheck["write_log_acts_on_after_step0_4"]
                                         and gcheck["graph_drug_nodes"] == gcheck["write_log_drugs_after_step0_4"])
    (OUT / "graph_check.json").write_text(json.dumps(gcheck, indent=2), encoding="utf-8")

    known = int(df.Action_direction.isin(["up", "down"]).sum())
    dirc = df.Action_direction.value_counts().reindex(["down", "up", "mixed", "unknown"], fill_value=0)
    dsrc = df.Direction_source.value_counts()
    numbers = {"seed_genes": len(seeds), "genes_with_drug": int(df.Gene_ID.nunique()), "unique_drugs": int(df.Drug_ID.nunique()),
               "links": len(df), "links_optimuskg": gcheck["edges_optimuskg"], "links_added": gcheck["edges_added"],
               "links_known_direction": known, "links_direction_comparable": int(df.direction_comparable.sum()),
               "relation_conflict_links": int(df.relation_conflict.astype("boolean").fillna(False).sum()),
               "approved_drugs": int(df.loc[df.is_approved.eq(True) | df.drug_stage.eq("approved"), "Drug_ID"].nunique())}
    pd.Series(numbers).to_csv(OUT / "summary_numbers.csv", header=["value"])
    by_stage = df.groupby("main_stage").agg(links=("Database_ID", "size"), genes=("Gene_ID", "nunique"),
                                           drugs=("Drug_ID", "nunique"),
                                           known_direction=("Action_direction", lambda s: int(s.isin(["up", "down"]).sum())),
                                           comparable=("direction_comparable", "sum"))
    by_stage.to_csv(OUT / "links_by_stage.csv")
    by_origin = df.groupby("Link_origin").agg(links=("Database_ID", "size"), genes=("Gene_ID", "nunique"),
                                             drugs=("Drug_ID", "nunique"),
                                             known_direction=("Action_direction", lambda s: int(s.isin(["up", "down"]).sum())))
    by_origin.to_csv(OUT / "links_by_origin.csv")
    raw = df.Action_raw.value_counts()

    def md(t):
        t = t.reset_index()
        return "\n".join(["| " + " | ".join(map(str, t.columns)) + " |", "|" + "---|" * t.shape[1]]
                         + ["| " + " | ".join(str(v) for v in r) + " |" for r in t.values])

    summary = f"""# Task 2: drug-gene links on the seed genes

Gene set: the {len(seeds)} seed genes (Task 1). Source: database paper2public, read-only
(`MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids`), queried {gcheck['queried']}.

- **{len(df):,} links** between **{df.Drug_ID.nunique():,} drugs** and **{df.Gene_ID.nunique()} of {len(seeds)} seed genes**
  ({gcheck['edges_optimuskg']:,} OptimusKG, {gcheck['edges_added']:,} added in step 0.4). No duplicate drug-gene pairs: {gcheck['duplicate_drug_gene_pairs'] == 0}.
- **Known drug direction: {known:,} links ({100 * known / len(df):.1f}%)**: down {dirc['down']:,}, up {dirc['up']:,};
  mixed {dirc['mixed']:,}, unknown {dirc['unknown']:,}.
- Direction source: {", ".join(f"{k} {v:,}" for k, v in dsrc.items())}.
- Links with a known drug direction AND a gene direction (comparable for the direction filter): {int(df.direction_comparable.sum()):,}.
  `0_genetic_risk` genes have no direction (`not_applicable`), so their links are never comparable.
- Graph totals: {gcheck['graph_drug_nodes']:,} Drug nodes, {gcheck['graph_acts_on_total']:,} ACTS_ON links; match the step 0.4 write log: {gcheck['graph_matches_write_log']}.

## By stage of the seed gene

{md(by_stage)}

## By link origin

{md(by_origin)}

## Most common actions

{md(raw.head(15).rename('links'))}

Direction rule: drug direction from the link's own action; OptimusKG actions INHIBITOR, ANTAGONIST, BLOCKER, NEGATIVE_(ALLOSTERIC_)MODULATOR,
ALLOSTERIC_ANTAGONIST, INVERSE_AGONIST, DEGRADER = down; AGONIST, PARTIAL_AGONIST, ACTIVATOR, POSITIVE_(ALLOSTERIC_)MODULATOR, OPENER,
STABILISER = up; MODULATOR, BINDING_AGENT, SUBSTRATE, RELEASING_AGENT = mixed; others unknown, then filled from step 0.5 where available.
Gene direction: increasing = positive, decreasing = negative.

Output: `drug_gene_edges.csv` ({len(df):,} rows), `summary_numbers.csv`, `links_by_stage.csv`, `links_by_origin.csv`, `graph_check.json`.
Script: `code/task2_drug_gene_edges.py`.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
