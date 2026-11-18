"""Task 4: drug coverage of the seed genes and repurposing candidate generation.

Main analysis uses pharmacological-target links only (same definition as Task 3). Any-link counts are kept as
secondary columns. Repurposing information (approval, phase, withdrawal, existing disease uses) comes from this
run's graph, plus the ChEMBL indications of the drugs added in step 0.4 (written as TREATS in step 0.8 and also
read from step 0.8's file, matched by name for EFO/MONDO ids outside the graph's liver set).

Candidate filter:
  phase >= 2  ->  not withdrawn  ->  no existing fibrosis-type liver use (graph TREATS/OFF_LABEL_FOR + ChEMBL
  indications)  ->  mechanism annotated (>= 1 link with a known direction or a stated mechanism)
Liver disease = EFO_0001421 'liver disease' + all PARENT_OF descendants, minus non-hepatic-tissue entries;
fibrosis-type = liver diseases matching fatty liver / steatohepatitis / cirrhosis / fibrosis / portal hypertension /
liver failure. Ordering (not a filter): histology seed genes hit, then target fraction, then all seed genes hit.
Gene annotation: each seed gene's most significant pathway from step 0.7 (top_enrichment_pathway).

Outputs: results/Task4/.  This is candidate generation, not ranking.

    python code/task4_drug_coverage.py
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, driver, now, read_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
EDGES = ROOT / "results/Task2/drug_gene_edges.csv"
TOP_PW = ROOT / "results/graph_build/step0_7_pathway_nafld/seed_gene_top_pathway.csv"
ADDED_IND = ROOT / "results/graph_build/step0_8_added_drug_indications/added_drug_indications_chembl.csv"
ADDED_IND_META = ROOT / "results/graph_build/step0_8_added_drug_indications/added_drug_indications_chembl.json"
OUT = ROOT / "results/Task4"
FIB_LIVER_NAME = (r"fatty liver|steatohepatitis|cirrhosis|(?:liver|hepatic) fibrosis|portal hypertension|"
                  r"(?:liver|hepatic) failure|NAFLD|NASH|MASLD|MASH")
PK_ACTIONS = ["ENZYME", "TRANSPORTER", "CARRIER"]
ADDED_LAYER = "added curated drug-target links"
LIVER_ROOT = "EFO_0001421"
LIVER_EXCLUDE = r"malaria|blackwater|gallbladder|gall bladder|ampulla|extrahepatic|bladder hepatoid|(?<!hepato)erythropoietic"
FIB_SPECTRUM = r"fatty liver|steatohepatitis|cirrhosis|fibrosis|portal hypertension|liver failure"
MIN_PHASE = 2
TOP_N_TABLE, TOP_N_NETWORK = 50, 30
STAGES = ["0_genetic_risk", "1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]
STAGE_COL = {"0_genetic_risk": "#7f8c8d", "1_steatosis": "#f1c40f", "2_NASH": "#e67e22", "3_fibrosis": "#c0392b",
             "4_HCC": "#6c3483", "drug": "#2e86c1"}


def md(df, cols):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join("" if pd.isna(r[c]) else str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    edges = pd.read_csv(EDGES, low_memory=False)
    edges["pharm"] = ((edges.Link_origin != "OptimusKG") | edges.external_action_direction.notna()
                      | (~edges.Action_raw.isin(PK_ACTIONS) & (edges.Action_raw != "unknown")))
    pe = edges[edges.pharm].copy()
    drug_ids = sorted(edges.Drug_ID.unique())

    drv = driver()
    with read_session(drv) as s:
        totals = pd.DataFrame(s.run("""UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:ACTS_ON]->(g:Gene)
            RETURN id AS Drug_ID, count(DISTINCT g) AS all_targets_graph,
                   count(DISTINCT CASE WHEN r.layer = $layer OR r.external_action_direction IS NOT NULL
                                            OR (r.action IS NOT NULL AND NOT r.action IN $pk) THEN g END) AS pharm_targets_graph""",
                                    ids=drug_ids, pk=PK_ACTIONS, layer=ADDED_LAYER).data())
        liver = pd.DataFrame(s.run("""MATCH (:Disease {node_id: $root})-[:PARENT_OF*0..20]->(c:Disease)
            RETURN DISTINCT c.node_id AS disease_id, c.name AS disease_name""", root=LIVER_ROOT).data())
        uses = pd.DataFrame(s.run("""UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:TREATS|OFF_LABEL_FOR|CONTRAINDICATED_IN]->(x:Disease)
            RETURN id AS Drug_ID, type(r) AS use_type, x.node_id AS disease_id, x.name AS disease_name,
                   coalesce(r.layer, 'OptimusKG') AS use_layer""", ids=drug_ids).data())
    drv.close()
    if liver.empty:
        raise SystemExit(f"no PARENT_OF tree under {LIVER_ROOT}")
    liver["disease_name"] = liver.disease_name.fillna(liver.disease_id)
    liver["excluded_non_hepatic"] = liver.disease_name.str.contains(LIVER_EXCLUDE, case=False, regex=True)
    liver["fibrosis_type"] = ~liver.excluded_non_hepatic & liver.disease_name.str.contains(FIB_SPECTRUM, case=False, regex=True)
    liver.sort_values(["excluded_non_hepatic", "disease_name"]).to_csv(OUT / "liver_disease_set.csv", index=False)
    liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
    fib_liver_ids = set(liver.loc[liver.fibrosis_type, "disease_id"])

    uses["use_source"] = np.where(uses.use_layer.eq("OptimusKG"), "graph (OptimusKG)", "graph (ChEMBL indication, step 0.8)")
    ind = pd.read_csv(ADDED_IND)
    ind_meta = json.loads(ADDED_IND_META.read_text(encoding="utf-8"))
    ind = ind[ind.Drug_ID.isin(drug_ids)]
    uses = pd.concat([uses, pd.DataFrame({
        "Drug_ID": ind.Drug_ID, "use_type": "TREATS", "disease_id": ind.efo_id_graph,
        "disease_name": ind.efo_term.fillna(ind.mesh_heading),
        "use_source": "ChEMBL indication file (phase " + ind.max_phase_for_ind.map(lambda x: f"{x:g}") + ")"})],
        ignore_index=True)
    uses["disease_name"] = uses.disease_name.fillna(uses.disease_id)
    chembl_use = uses.use_source.str.contains("ChEMBL")
    name_fib = uses.disease_name.str.contains(FIB_LIVER_NAME, case=False, regex=True)
    uses["is_liver"] = uses.disease_id.isin(liver_ids) | (chembl_use & name_fib)
    uses["is_fib_liver"] = uses.disease_id.isin(fib_liver_ids) | (chembl_use & name_fib)

    tp = pd.read_csv(TOP_PW)
    top_pw = dict(zip(tp.ensembl_id, (tp.top_pathway_db + ": " + tp.top_pathway).where(tp.top_pathway.notna())))
    seeds["top_enrichment_pathway"] = seeds.ensembl_id.map(top_pw)

    def per_drug(e):
        g = e.groupby("Drug_ID")
        out = pd.DataFrame({"n_seed_genes": g.Gene_ID.nunique(),
                            "n_histology_genes": e[e.main_stage != "0_genetic_risk"].groupby("Drug_ID").Gene_ID.nunique()})
        for st in STAGES:
            out[st] = e[e.main_stage == st].groupby("Drug_ID").Gene_ID.nunique()
        for st in STAGES[1:]:
            for d, lab in [("increasing", "up"), ("decreasing", "down")]:
                out[f"{st}_{lab}"] = e[(e.main_stage == st) & (e.seed_direction == d)].groupby("Drug_ID").Gene_ID.nunique()
        for d in ["down", "up", "mixed", "unknown"]:
            out[f"n_{d}_actions"] = e[e.Action_direction == d].groupby("Drug_ID").size()
        out["mechanism_annotated"] = g.apply(lambda s: bool((s.Action_direction != "unknown").any()
                                                            or s.mechanism_of_action.notna().any()), include_groups=False)
        out["target_genes"] = g.Gene.agg(lambda s: ";".join(sorted(set(s))))
        return out

    main_ = per_drug(pe)
    anyl = edges.groupby("Drug_ID").Gene_ID.nunique().rename("n_seed_genes_any_link")
    info = (edges.sort_values("Drug_ID").drop_duplicates("Drug_ID")
            .set_index("Drug_ID")[["Drug", "molecule_type", "max_phase", "is_approved", "withdrawn", "drug_stage"]])
    info["added_drug_node"] = info.is_approved.isna() & info.drug_stage.notna()
    info["is_approved"] = info.is_approved.astype("boolean").fillna(info.drug_stage.eq("approved")).astype(bool)
    info["phase_from_stage"] = info.added_drug_node & info.max_phase.isna() & info.drug_stage.eq("approved")
    info.loc[info.phase_from_stage, "max_phase"] = 4.0
    info = info.drop(columns="drug_stage")

    cov = info.join(main_, how="left").join(anyl).join(totals.set_index("Drug_ID"))
    txt = {"Drug", "molecule_type", "max_phase", "is_approved", "withdrawn", "added_drug_node", "phase_from_stage",
           "target_genes", "mechanism_annotated"}
    num = [c for c in cov.columns if c not in txt]
    cov["mechanism_annotated"] = cov.mechanism_annotated.astype("boolean").fillna(False).astype(bool)
    cov[num] = cov[num].fillna(0).astype(int)
    cov["seed_target_fraction"] = np.where(cov.pharm_targets_graph > 0,
                                           cov.n_seed_genes / cov.pharm_targets_graph.replace(0, np.nan), np.nan).round(3)
    use_ok = uses.use_type.isin(["TREATS", "OFF_LABEL_FOR"])
    cov["n_treats"] = uses[uses.use_type == "TREATS"].groupby("Drug_ID").disease_id.nunique()
    cov["has_liver_indication"] = uses[uses.is_liver & use_ok].groupby("Drug_ID").size() > 0
    cov["liver_indications"] = uses[uses.is_liver & use_ok].groupby("Drug_ID").disease_name.agg(lambda s: "; ".join(sorted(set(s))))
    cov["has_fibrosis_type_liver_use"] = uses[uses.is_fib_liver & use_ok].groupby("Drug_ID").size() > 0
    cov["fibrosis_type_liver_uses"] = uses[uses.is_fib_liver & use_ok].groupby("Drug_ID").disease_name.agg(lambda s: "; ".join(sorted(set(s))))
    ci = uses.is_liver & (uses.use_type == "CONTRAINDICATED_IN")
    cov["liver_contraindicated"] = uses[ci].groupby("Drug_ID").size() > 0
    cov["liver_contraindications"] = uses[ci].groupby("Drug_ID").disease_name.agg(lambda s: "; ".join(sorted(set(s))))
    cov["current_indications_examples"] = uses[uses.use_type == "TREATS"].groupby("Drug_ID").disease_name.agg(
        lambda s: "; ".join(sorted(set(s))[:5]))
    cov["n_treats"] = cov.n_treats.fillna(0).astype(int)
    for c in ["has_liver_indication", "has_fibrosis_type_liver_use", "liver_contraindicated"]:
        cov[c] = cov[c].astype("boolean").fillna(False).astype(bool)
    cov["kept_with_other_liver_use"] = cov.has_liver_indication & ~cov.has_fibrosis_type_liver_use
    cov = cov.reset_index().rename(columns={"index": "Drug_ID"})
    stage_cols = STAGES + [f"{st}_{lab}" for st in STAGES[1:] for lab in ["up", "down"]]
    order = (["Drug", "Drug_ID", "molecule_type", "max_phase", "is_approved", "withdrawn", "added_drug_node", "phase_from_stage",
              "n_seed_genes", "n_histology_genes"] + stage_cols
             + ["mechanism_annotated", "n_down_actions", "n_up_actions", "n_mixed_actions", "n_unknown_actions",
                "pharm_targets_graph", "all_targets_graph", "seed_target_fraction", "n_seed_genes_any_link", "n_treats",
                "has_liver_indication", "liver_indications", "has_fibrosis_type_liver_use", "fibrosis_type_liver_uses",
                "kept_with_other_liver_use", "liver_contraindicated", "liver_contraindications",
                "current_indications_examples", "target_genes"])
    cov = cov[order].sort_values(["n_seed_genes", "n_histology_genes", "seed_target_fraction"], ascending=False)
    cov.to_csv(OUT / "drug_coverage.csv", index=False)
    cov[cov.n_seed_genes > 0].head(TOP_N_TABLE).to_csv(OUT / "top50_by_coverage.csv", index=False)

    # ---------------------------------------------------------- candidate funnel
    hit = cov[cov.n_seed_genes >= 1]
    f1 = hit[hit.max_phase >= MIN_PHASE]
    f2 = f1[~f1.withdrawn.astype("boolean").fillna(False).astype(bool)]
    f3 = f2[~f2.has_fibrosis_type_liver_use]
    f4 = f3[f3.mechanism_annotated]
    funnel = pd.DataFrame([
        ("Drugs with >= 1 link to a seed gene (any action)", cov.Drug_ID.nunique()),
        ("... with >= 1 pharmacological-target link (coverage table)", len(hit)),
        (f"... max phase >= {MIN_PHASE}", len(f1)),
        ("... and not withdrawn", len(f2)),
        ("... and no existing fibrosis-type liver use = repurposing candidates", len(f3)),
        ("... and mechanism annotated = main list", len(f4))], columns=["step", "drugs"])
    funnel["removed_at_step"] = [0, funnel.drugs[0] - funnel.drugs[1]] + list(-funnel.drugs.diff().iloc[2:].astype(int))
    funnel.to_csv(OUT / "candidate_funnel.csv", index=False)
    cand = f3.sort_values(["n_histology_genes", "seed_target_fraction", "n_seed_genes"], ascending=False).copy()
    cand.insert(0, "candidate_order", range(1, len(cand) + 1))
    cand.to_csv(OUT / "repurposing_candidates.csv", index=False)
    mech = cand[cand.mechanism_annotated].copy()
    mech.insert(1, "mechanism_order", range(1, len(mech) + 1))
    mech.to_csv(OUT / "repurposing_candidates_mechanism_annotated.csv", index=False)

    # ---------------------------------------------------------- network (top 30 mechanism-annotated candidates)
    top = mech.head(TOP_N_NETWORK)
    ne = pe[pe.Drug_ID.isin(top.Drug_ID)]
    gann = seeds.set_index("ensembl_id")
    pd.DataFrame({"source": ne.Drug_ID, "target": ne.Gene_ID, "source_name": ne.Drug, "target_name": ne.Gene,
                  "action": ne.Action_raw, "action_direction": ne.Action_direction, "evidence_level": ne.evidence_level,
                  "link_origin": ne.Link_origin, "gene_stage": ne.main_stage, "gene_direction": ne.gene_direction}).to_csv(
        OUT / "network_edges.csv", index=False)
    gn = gann.loc[ne.Gene_ID.unique()]
    nodes = pd.concat([pd.DataFrame({"id": top.Drug_ID, "label": top.Drug, "node_type": "drug", "stage": "drug",
                                     "direction": np.nan, "top_enrichment_pathway": np.nan, "is_approved": top.is_approved,
                                     "n_seed_genes": top.n_seed_genes}),
                       pd.DataFrame({"id": gn.index, "label": gn.gene, "node_type": "gene", "stage": gn.main_stage,
                                     "direction": gn.direction, "top_enrichment_pathway": gn.top_enrichment_pathway,
                                     "is_approved": np.nan, "n_seed_genes": np.nan})], ignore_index=True)
    nodes.to_csv(OUT / "network_nodes.csv", index=False)
    G = nx.Graph()
    for _, r in nodes.iterrows():
        G.add_node(r.id, label=r.label, node_type=r.node_type, stage=r.stage)
    for _, r in ne.iterrows():
        G.add_edge(r.Drug_ID, r.Gene_ID, direction=r.Action_direction)
    if len(G):
        pos = nx.spring_layout(G, k=0.9, seed=42, iterations=300)
        fig, ax = plt.subplots(figsize=(15, 12))
        for d, c in {"down": "#2471a3", "up": "#c0392b", "mixed": "#d4ac0d", "unknown": "#c5ccd3"}.items():
            es = [(u, v) for u, v, a in G.edges(data=True) if a["direction"] == d]
            nx.draw_networkx_edges(G, pos, edgelist=es, edge_color=c, width=1.6 if d in ("down", "up") else 0.6, ax=ax)
        gnodes = [n for n, a in G.nodes(data=True) if a["node_type"] == "gene"]
        dnodes = [n for n, a in G.nodes(data=True) if a["node_type"] == "drug"]
        nx.draw_networkx_nodes(G, pos, nodelist=gnodes, node_size=120, node_color=[STAGE_COL[G.nodes[n]["stage"]] for n in gnodes], ax=ax)
        nx.draw_networkx_nodes(G, pos, nodelist=dnodes, node_shape="s", node_size=380, node_color=STAGE_COL["drug"], ax=ax)
        nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]["label"] for n in G.nodes}, font_size=7, ax=ax)
        handles = [plt.Line2D([], [], marker="s", ls="", color=STAGE_COL["drug"], ms=10, label="Drug (candidate)")]
        handles += [plt.Line2D([], [], marker="o", ls="", color=STAGE_COL[s_], ms=8, label=f"Seed gene, {s_}") for s_ in STAGES]
        handles += [plt.Line2D([], [], color=c, lw=2, label=l) for c, l in
                    [("#2471a3", "drug pushes gene down"), ("#c0392b", "drug pushes gene up"), ("#c5ccd3", "direction unknown")]]
        ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=False)
        ax.set_title(f"Top {TOP_N_NETWORK} mechanism-annotated repurposing candidates and their seed-gene targets",
                     loc="left", fontweight="bold")
        ax.axis("off")
        fig.savefig(OUT / "network_top30_candidates.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    added_all = cov[cov.added_drug_node & (cov.n_seed_genes > 0)]
    summary = f"""# Task 4: drug coverage and repurposing candidates

Seed genes: 280 (Task 1). Drug links: Task 2 (pharmacological-target links for the main analysis). Generated {now()}.

## Candidate funnel

{md(funnel, ['step', 'drugs', 'removed_at_step'])}

- Liver diseases: `{LIVER_ROOT}` + PARENT_OF descendants = {len(liver)} terms, minus {int(liver.excluded_non_hepatic.sum())}
  non-hepatic-tissue entries = {len(liver_ids)} liver diseases; {len(fib_liver_ids)} are fibrosis-type (`liver_disease_set.csv`).
- Drugs hitting seed genes with any liver TREATS/OFF_LABEL_FOR use: {int(hit.has_liver_indication.sum())};
  with a fibrosis-type liver use (excluded): {int(hit.has_fibrosis_type_liver_use.sum())};
  kept although they have another liver use: {int(f3.kept_with_other_liver_use.sum())}.
- Liver contraindication (kept here, removed in the ranking step): {int(f4.liver_contraindicated.sum())} of the {len(f4)} final candidates ({int(f3.liver_contraindicated.sum())} of the {len(f3)} drugs before the mechanism filter).
- Drugs added in step 0.4 among the coverage table: {len(added_all)}; approved with no phase number (treated as phase 4):
  {int(added_all.phase_from_stage.sum())}; no phase number at all (cannot confirm phase >= 2, excluded): {int(added_all.max_phase.isna().sum())}.
  Their existing uses come from ChEMBL indications ({ind_meta['chembl_release']}, {ind_meta['chembl_release_date']}, retrieved {ind_meta['retrieved']}).

## Main list: first 20 mechanism-annotated candidates (order = histology seed genes hit, then target fraction)

{md(mech.head(20), ['mechanism_order', 'Drug', 'max_phase', 'n_seed_genes', 'n_histology_genes', '3_fibrosis', 'n_down_actions', 'n_up_actions', 'target_genes'])}

This is candidate generation, not a ranking of drug quality.

Outputs: `drug_coverage.csv` (every drug with a seed-gene link), `top50_by_coverage.csv`, `candidate_funnel.csv`,
`repurposing_candidates.csv`, `repurposing_candidates_mechanism_annotated.csv` (main list), `liver_disease_set.csv`,
`network_nodes.csv` / `network_edges.csv` / `network_top30_candidates.png`. Gene annotation: `top_enrichment_pathway`
(most significant step 0.7 pathway per seed gene). Script: `code/task4_drug_coverage.py`.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
