"""Independent pharmacological and mechanistic evaluation, steps 1 and 2.

Step 1  export the 116 conventional-MPC drugs, identify the 38 retained by the direction filter, verify the 78 excluded
        (reason for each) and mark the annotation-only group (decidability-only filter).
Step 2  extract every drug-target pair (all ACTS_ON targets, not only seed targets) of the 116 drugs with the
        pharmacological action label, the derived direction, sources, evidence level and the seed-gene status.

Reads results of the main pipeline and the graph (read-only Neo4j). Writes results/Evaluation/mechanistic/.

    NEO4J_PASSWORD=... NEO4J_DATABASE=neo4j python code/mech_step1_2_export.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, driver, read_session  # noqa: E402
from task2_drug_gene_edges import edge_direction  # noqa: E402
from evalcore import Data, Routes, adjacency, path_scores, pathway_links, start_links  # noqa: E402

OUT = ROOT / "results/Evaluation/mechanistic"
OUT.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ step 1
base = pd.read_csv(ROOT / "results/Baseline_no_direction/RunA_all_seed_genes_FINAL_ranked_drugs.csv")
main = pd.read_csv(ROOT / "results/Ranking/RunA_all_seed_genes_FINAL_ranked_drugs.csv")
summ = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_drugs_direction_summary.csv").set_index("Drug_ID")

D = Data()
A = adjacency(D)
linked = pathway_links(D, D.seeds)
R = Routes(D, A, D.seeds.i.to_numpy(), linked)
mpc = path_scores(D, R, start_links(D, D.seeds.i.to_numpy(), "none"))
ann = path_scores(D, R, start_links(D, D.seeds.i.to_numpy(), "decidable"))
prop = path_scores(D, R, start_links(D, D.seeds.i.to_numpy(), "proposed"))

ids_mpc = set(mpc[mpc > 0].index)
assert ids_mpc == set(base.Drug_ID), "offline MPC list differs from the Neo4j baseline list"
assert set(prop[prop > 0].index) == set(main.Drug_ID), "offline proposed list differs from the Neo4j main list"
retained = set(main.Drug_ID)
assert retained <= ids_mpc
print("MPC drugs", len(ids_mpc), "| retained", len(retained), "| excluded", len(ids_mpc - retained))

t = base[["Rank", "Drug", "Drug_ID", "BB routes to NAFLD", "B routes"]].rename(columns={"Rank": "mpc_rank", "BB routes to NAFLD": "mpc_BB_routes", "B routes": "mpc_B_routes"})
t["direction_group"] = t.Drug_ID.map(summ.direction_group)
t["n_consistent_links"] = t.Drug_ID.map(summ.n_consistent)
t["n_inconsistent_links"] = t.Drug_ID.map(summ.n_inconsistent)
t["n_undecided_links"] = t.Drug_ID.map(summ.n_not_decided)
t["retained_by_filter"] = t.Drug_ID.isin(retained)
t["proposed_rank"] = t.Drug_ID.map(main.set_index("Drug_ID").Rank)
t["proposed_BB_routes"] = t.Drug_ID.map(main.set_index("Drug_ID")["BB routes to NAFLD"])
t["annotation_only_group"] = t.Drug_ID.isin(set(ann[ann > 0].index))      # decidability-only filter, drugs with routes
t["liver_contraindicated"] = t.Drug_ID.map(summ.liver_contraindicated)


def reason(r):
    if r.retained_by_filter:
        return "retained"
    if r.direction_group == "no decided link":
        return "excluded: no decided drug direction on any seed link"
    if r.direction_group == "only inconsistent":
        return "excluded: all decided links inconsistent with the disease direction"
    if r.direction_group == "only consistent":
        return "excluded: direction-consistent, but its consistent targets have no B route, so BB is never searched (its MPC B routes come from undecided targets)"
    return "excluded: " + str(r.direction_group)


t["group_reason"] = t.apply(reason, axis=1)
t["group"] = np.where(t.retained_by_filter, "retained_38", "excluded_78")
t = t.sort_values(["retained_by_filter", "mpc_rank", "Drug"], ascending=[False, True, True])
t.to_csv(OUT / "step1_drug_groups.csv", index=False)
chk = t.group_reason.value_counts()
print(chk.to_string())
assert (t.group == "retained_38").sum() == 38 and (t.group == "excluded_78").sum() == 78

# ------------------------------------------------------------------ step 2
QUERY = """
UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:ACTS_ON]->(g:Gene)
RETURN d.node_id AS Drug_ID, coalesce(d.recovered_name, d.name) AS Drug, g.node_id AS Gene_ID, g.symbol AS Gene,
       r.layer AS layer, r.action AS okg_action, r.relation_assertions AS relation_assertions, r.relation_conflict AS relation_conflict,
       r.sources_direct AS sources_direct, r.evidence_level AS evidence_level, r.mechanisms_of_action AS okg_mechanism,
       r.action_labels AS ext_action_labels, r.push AS ext_push, r.sources AS ext_sources, r.evidence AS ext_evidence,
       r.mechanism AS ext_mechanism, r.external_action_direction AS external_action_direction,
       r.external_action_label AS external_action_label, r.external_action_sources AS external_action_sources
"""
drv = driver()
with read_session(drv) as s:
    p = pd.DataFrame(s.run(QUERY, ids=sorted(ids_mpc)).data())
drv.close()


def j(v):
    return ";".join(map(str, v)) if isinstance(v, list) else v


def blank(x):
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == ""


for c in ["sources_direct", "okg_mechanism", "relation_assertions", "ext_sources"]:
    p[c] = p[c].map(j)
added = p.layer.eq("added curated drug-target links")
p["link_origin"] = np.where(added, "added (ChEMBL/DGIdb/GtoPdb/OpenTargets)", "OptimusKG")
p["action_raw"] = [(lab if not blank(lab) else "unknown") if a else (act if not blank(act) else "unknown")
                   for a, lab, act in zip(added, p.ext_action_labels, p.okg_action)]
dirn = [({"down": "down", "up": "up", "conflicting": "mixed"}.get(str(pu), "unknown")) if a else edge_direction(act, ass)
        for a, pu, act, ass in zip(added, p.ext_push, p.okg_action, p.relation_assertions)]
p["action_direction"] = dirn
fill = (~added) & p.action_direction.eq("unknown") & p.external_action_direction.isin(["down", "up"])
p.loc[fill, "action_direction"] = p.loc[fill, "external_action_direction"]
p["direction_source"] = np.select([added & p.action_direction.ne("unknown"), fill, p.action_direction.ne("unknown")],
                                  ["added link", "external look-up (GtoPdb/DGIdb/ChEMBL)", "OptimusKG action"], "none")
p["mechanism_text"] = p.okg_mechanism.where(~added, p.ext_mechanism)
p["pair_sources"] = np.where(added, p.ext_sources.astype(str), p.sources_direct.astype(str))
sd = D.seeds.set_index("ensembl_id")
p["is_seed_gene"] = p.Gene_ID.isin(sd.index)
p["seed_direction"] = p.Gene_ID.map(sd.direction)
p["seed_stage"] = p.Gene_ID.map(sd.main_stage)
a, g_ = p.action_direction, p.seed_direction
p["direction_check"] = np.select([p.is_seed_gene & (((a == "down") & (g_ == "increasing")) | ((a == "up") & (g_ == "decreasing"))),
                                  p.is_seed_gene & (((a == "down") & (g_ == "decreasing")) | ((a == "up") & (g_ == "increasing"))),
                                  p.is_seed_gene], ["consistent", "inconsistent", "undecided"], "not a seed gene")
p["group"] = p.Drug_ID.map(t.set_index("Drug_ID").group)
p["mpc_rank"] = p.Drug_ID.map(t.set_index("Drug_ID").mpc_rank)
cols = ["group", "mpc_rank", "Drug_ID", "Drug", "Gene_ID", "Gene", "is_seed_gene", "seed_stage", "seed_direction", "action_raw",
        "action_direction", "direction_source", "direction_check", "link_origin", "evidence_level", "relation_conflict",
        "pair_sources", "mechanism_text", "external_action_label", "external_action_sources"]
p = p[cols].sort_values(["group", "mpc_rank", "Drug", "is_seed_gene", "Gene"], ascending=[False, True, True, False, True])
p.to_csv(OUT / "step2_drug_target_pairs.csv", index=False)

n_pairs = len(p)
sp = p[p.is_seed_gene]
summary = pd.DataFrame([
    {"group": g, "drugs": x.Drug_ID.nunique(), "pairs": len(x), "seed_pairs": int(x.is_seed_gene.sum()),
     "pairs_with_direction": int(x.action_direction.isin(["down", "up"]).sum()),
     "distinct_targets": x.Gene_ID.nunique(), "distinct_seed_targets": x[x.is_seed_gene].Gene_ID.nunique()}
    for g, x in p.groupby("group")])
summary.to_csv(OUT / "step2_summary.csv", index=False)
print(summary.to_string(index=False))
print("pairs total", n_pairs, "| seed pairs", len(sp), "| drugs without any pair:", len(ids_mpc - set(p.Drug_ID)))
print(p.action_raw.value_counts().head(12).to_string())
