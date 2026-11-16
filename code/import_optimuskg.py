"""Step 0.0: build the graph fresh from the OptimusKG parquet files (Dataverse doi:10.7910/DVN/IYNGEV).

Nodes (label <- parquet)               key property node_id = OptimusKG id
  Gene <- nodes/gene                   ensembl_id, symbol, name, biotype
  Drug <- nodes/drug                   name, molecule_type, max_phase, is_approved, has_been_withdrawn,
                                       smiles, inchi_key, cd_formula, cd_mol_weight, year_of_first_approval,
                                       black_box_warning, sources_direct, source_ids (+ *_source = 'OptimusKG')
  Disease <- nodes/disease             name, code
  Pathway <- nodes/pathway             name, source = 'Reactome'
  BiologicalProcess <- nodes/biological_process   name
  Phenotype <- nodes/phenotype         name

Relationships (type <- parquet, relation)
  (Drug)-[:ACTS_ON]->(Gene)            drug_gene; action = relation, plus sources_direct/indirect,
                                       relation_assertions, relation_conflict, source_ids, source_urls,
                                       mechanisms_of_action, n_sources, evidence_level
  (Gene)-[:INTERACTS_WITH]->(Gene)     gene_gene (undirected in OptimusKG; stored once, queried undirected)
  (Gene)-[:IN_PATHWAY]->(Pathway)      pathway_gene
  (Drug)-[:TREATS]->(Disease)          drug_disease INDICATION
  (Drug)-[:OFF_LABEL_FOR]->(Disease)   drug_disease OFF_LABEL_USE
  (Drug)-[:CONTRAINDICATED_IN]->(Disease)  drug_disease CONTRAINDICATION
  (Gene)-[:ASSOCIATED_WITH]->(Disease) disease_gene; evidence_score, evidence_count
  (Disease)-[:PARENT_OF]->(Disease)    disease_disease PARENT (from = parent)
  (Disease)-[:HAS_PHENOTYPE]->(Phenotype)  disease_phenotype PHENOTYPE_PRESENT
  (Gene)-[:PARTICIPATES_IN]->(BiologicalProcess)  biological_process_gene

Not imported: anatomy (incl. the BGEE expression layer anatomy_gene), cellular component,
molecular function, exposure, drug_drug, drug_phenotype, drug_biological_process, phenotype_gene
and the ontology-internal edges other than disease_disease. No expression or progression layer.

Dry run by default (reads parquet, checks every edge endpoint exists, prints the plan).
--apply requires the target database to be empty, writes, then checks graph counts equal the plan.

    python code/import_optimuskg.py            # dry run
    python code/import_optimuskg.py --apply
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, TARGET_DB, counts, driver, now, run_batches, write_log, write_session  # noqa: E402

OKG = ROOT / "data/optimuskg"
OUT = ROOT / "results/graph_build/step0_0_import"


def norm(v):
    """parquet value -> Neo4j-storable value (None, scalar, or list of scalars)."""
    if v is None:
        return None
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if isinstance(v, (list, tuple)):
        out = [norm(x) for x in v]
        out = [x for x in out if x is not None and not isinstance(x, (dict, list))]
        return out or None
    if isinstance(v, dict):
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def src(p, which):
    s = (p or {}).get("sources") or {}
    return norm(s.get(which))


def assertions(p):
    a = (p or {}).get("relation_assertions")
    if a is None:
        return None
    return ";".join(f"{x.get('source')}:{x.get('relation')}" for x in a) or None


def clean(d):
    return {k: v for k, v in d.items() if v is not None}


# ---------------------------------------------------------------- nodes
def load_nodes():
    nodes = {}

    g = pd.read_parquet(OKG / "nodes/gene.parquet")
    nodes["Gene"] = [clean({"node_id": i, "ensembl_id": i, "symbol": norm(p.get("symbol")),
                            "name": norm(p.get("name")), "biotype": norm(p.get("biotype"))})
                     for i, p in zip(g["id"], g["properties"])]

    d = pd.read_parquet(OKG / "nodes/drug.parquet")
    drugs = []
    for i, p in zip(d["id"], d["properties"]):
        r = {"node_id": i, "name": norm(p.get("name")), "molecule_type": norm(p.get("type")),
             "max_phase": norm(p.get("maximum_clinical_trial_phase")),
             "is_approved": norm(p.get("is_approved")), "has_been_withdrawn": norm(p.get("has_been_withdrawn")),
             "smiles": norm(p.get("canonical_smiles")), "inchi_key": norm(p.get("inchi_key")),
             "cd_formula": norm(p.get("cd_formula")), "cd_mol_weight": norm(p.get("cd_mol_weight")),
             "year_of_first_approval": norm(p.get("year_of_first_approval")),
             "black_box_warning": norm(p.get("black_box_warning")),
             "sources_direct": src(p, "direct"), "source_ids": norm(p.get("source_ids"))}
        for f in ("molecule_type", "max_phase", "smiles"):
            if r[f] is not None:
                r[f + "_source"] = "OptimusKG"
        # OptimusKG stores an unnamed drug's name as its accession: that is "no name"
        r["named_in_optimuskg"] = bool(r["name"]) and r["name"] != i
        drugs.append(clean(r))
    nodes["Drug"] = drugs

    for label, f, extra in [("Disease", "disease", ["code"]), ("Pathway", "pathway", []),
                            ("BiologicalProcess", "biological_process", []), ("Phenotype", "phenotype", [])]:
        t = pd.read_parquet(OKG / f"nodes/{f}.parquet")
        rows = []
        for i, p in zip(t["id"], t["properties"]):
            r = {"node_id": i, "name": norm(p.get("name"))}
            for e in extra:
                r[e] = norm(p.get(e))
            if label == "Pathway":
                r["source"] = "Reactome"
            rows.append(clean(r))
        nodes[label] = rows
    return nodes


# ---------------------------------------------------------------- edges
# file, relation filter, rel type, start label, end label, reverse (parquet 'to' is our start)
EDGE_SPECS = [
    ("drug_gene", None, "ACTS_ON", "Drug", "Gene", False),
    ("gene_gene", None, "INTERACTS_WITH", "Gene", "Gene", False),
    ("pathway_gene", None, "IN_PATHWAY", "Gene", "Pathway", True),
    ("drug_disease", "INDICATION", "TREATS", "Drug", "Disease", False),
    ("drug_disease", "OFF_LABEL_USE", "OFF_LABEL_FOR", "Drug", "Disease", False),
    ("drug_disease", "CONTRAINDICATION", "CONTRAINDICATED_IN", "Drug", "Disease", False),
    ("disease_gene", None, "ASSOCIATED_WITH", "Gene", "Disease", True),
    ("disease_disease", None, "PARENT_OF", "Disease", "Disease", False),
    ("disease_phenotype", None, "HAS_PHENOTYPE", "Disease", "Phenotype", False),
    ("biological_process_gene", None, "PARTICIPATES_IN", "Gene", "BiologicalProcess", True),
]


def edge_props(fname, p):
    p = p or {}
    if fname == "drug_gene":
        a = assertions(p)
        n = len([x for x in (a or "").split(";") if x])
        conflict = bool(p.get("relation_conflict"))
        return clean({"sources_direct": src(p, "direct"), "sources_indirect": src(p, "indirect"),
                      "relation_assertions": a, "relation_conflict": conflict,
                      "source_ids": norm(p.get("source_ids")), "source_urls": norm(p.get("source_urls")),
                      "mechanisms_of_action": norm(p.get("mechanisms_of_action")), "n_sources": n,
                      "evidence_level": ("conflicting" if conflict else "multi-source" if n > 1 else
                                         "single-source" if n == 1 else "unattributed")})
    if fname == "drug_disease":
        return clean({"sources_direct": src(p, "direct"), "sources_indirect": src(p, "indirect"),
                      "relation_assertions": assertions(p), "relation_conflict": bool(p.get("relation_conflict")),
                      "highest_clinical_trial_phase": norm(p.get("highest_clinical_trial_phase")),
                      "reference_ids": norm(p.get("reference_ids"))})
    if fname == "disease_gene":
        return clean({"evidence_score": norm(p.get("evidence_score")),
                      "evidence_count": norm(p.get("evidence_count")), "sources_direct": src(p, "direct")})
    return clean({"sources_direct": src(p, "direct"), "sources_indirect": src(p, "indirect")})


def load_edges(node_ids):
    """Select rows per relationship type and check endpoints. Returns {type: (spec, DataFrame)} and the plan."""
    plan, edges, cache = [], {}, {}
    for spec in EDGE_SPECS:
        fname, rel_filter, rtype, a, b, rev = spec
        if fname not in cache:
            cache = {fname: pd.read_parquet(OKG / f"edges/{fname}.parquet")}
        t = cache[fname]
        sub = t if rel_filter is None else t[t.relation == rel_filter]
        s, e = (sub["to"], sub["from"]) if rev else (sub["from"], sub["to"])
        ok = s.isin(node_ids[a]).to_numpy() & e.isin(node_ids[b]).to_numpy()
        keep = pd.DataFrame({"s": s[ok].to_numpy(), "e": e[ok].to_numpy(),
                             "relation": sub["relation"][ok].to_numpy(), "properties": sub["properties"][ok].to_numpy()})
        edges[rtype] = (spec, keep)
        plan.append({"file": f"edges/{fname}.parquet", "relation": rel_filter or "all", "file_rows": len(t),
                     "selected_rows": len(sub), "missing_endpoint_dropped": int((~ok).sum()),
                     "type": rtype, "pattern": f"({a})-[:{rtype}]->({b})", "planned": len(keep)})
        print(f"  {rtype:20s} {len(keep):>10,} of {len(sub):>10,} rows ({int((~ok).sum())} with a missing endpoint)",
              flush=True)
    return edges, plan


def edge_rows(fname, df):
    for si, ei, rel, p in zip(df["s"], df["e"], df["relation"], df["properties"]):
        props = edge_props(fname, p)
        if fname == "drug_gene":
            props["action"] = rel
        yield {"s": si, "e": ei, "p": props}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    print("reading nodes ...", flush=True)
    nodes = load_nodes()
    node_ids = {k: set(r["node_id"] for r in v) for k, v in nodes.items()}
    for k, v in nodes.items():
        dup = len(v) - len(node_ids[k])
        print(f"  {k:18s} {len(v):>8,}" + (f"  DUPLICATE IDS: {dup}" if dup else ""))
        if dup:
            raise SystemExit("duplicate node ids")
    print("reading edges ...", flush=True)
    edges, plan = load_edges(node_ids)
    expected = {"nodes": {k: len(v) for k, v in nodes.items()}, "relationships": {k: len(v[1]) for k, v in edges.items()}}
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(plan).to_csv(OUT / "edge_plan.csv", index=False)
    (OUT / "expected_counts.json").write_text(json.dumps(expected, indent=2), encoding="utf-8")
    if not args.apply:
        print("dry run: nothing written. Plan in results/graph_build/step0_0_import/")
        return

    drv = driver()
    with write_session(drv) as s:
        before = counts(s)
        if before["nodes_total"] or before["relationships_total"]:
            raise SystemExit(f"{TARGET_DB} is not empty ({before['nodes_total']} nodes); import needs an empty database.")
        started = now()
        for label in nodes:
            s.run(f"CREATE CONSTRAINT {label.lower()}_node_id IF NOT EXISTS FOR (n:{label}) REQUIRE n.node_id IS UNIQUE")
        s.run("CREATE INDEX gene_symbol IF NOT EXISTS FOR (n:Gene) ON (n.symbol)")
        s.run("CALL db.awaitIndexes(600)")
        for label, rows in nodes.items():
            n = run_batches(s, f"UNWIND $rows AS r CREATE (n:{label}) SET n = r RETURN count(n)", rows, batch=5000)
            print(f"  wrote {label}: {n:,}", flush=True)
        for rtype, ((fname, _, _, a, b, _), df) in edges.items():
            q = (f"UNWIND $rows AS r MATCH (a:{a} {{node_id: r.s}}) MATCH (b:{b} {{node_id: r.e}}) "
                 f"CREATE (a)-[x:{rtype}]->(b) SET x = r.p RETURN count(x)")
            n, buf = 0, []
            for row in edge_rows(fname, df):
                buf.append(row)
                if len(buf) == 20000:
                    n += run_batches(s, q, buf, batch=20000)
                    buf = []
            if buf:
                n += run_batches(s, q, buf, batch=20000)
            print(f"  wrote {rtype}: {n:,}", flush=True)
        after = counts(s)
    drv.close()

    mism = {k: (v, after["nodes"].get(k)) for k, v in expected["nodes"].items() if after["nodes"].get(k) != v}
    mism.update({k: (v, after["relationships"].get(k)) for k, v in expected["relationships"].items()
                 if after["relationships"].get(k) != v})
    extra = set(after["relationships"]) - set(expected["relationships"])
    (OUT / "undo_import.cypher").write_text(
        "// Undo step 0.0: the database was empty before the import, so undo = delete everything.\n"
        f"// Run only against {TARGET_DB}.\n"
        "MATCH (n) CALL (n) { DETACH DELETE n } IN TRANSACTIONS OF 10000 ROWS;\n", encoding="utf-8")
    write_log(OUT / "write_log.json", {"step": "0.0 import OptimusKG", "database": TARGET_DB, "started": started,
                                       "finished": now(), "before": before, "after": after, "expected": expected,
                                       "mismatches": mism, "unexpected_types": sorted(extra),
                                       "undo": "undo_import.cypher"})
    print(json.dumps(after, indent=1))
    if mism or extra:
        raise SystemExit(f"COUNT CHECK FAILED: {mism} {extra}")
    print("count check passed: graph counts equal the plan")


if __name__ == "__main__":
    main()
