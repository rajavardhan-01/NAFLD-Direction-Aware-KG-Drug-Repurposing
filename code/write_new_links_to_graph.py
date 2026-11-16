"""Step 0.4 (write): add the curated new drug-gene links as Drug -ACTS_ON-> Gene.

Input: results/graph_build/step0_4_new_links/new_links_medicines_curated.csv (from new_links_search.py).
Rules: one relationship per curated drug-gene pair, CREATE (never MERGE),
existing evidence untouched. A drug key is matched to an existing Drug node by its own ChEMBL ID or any
ChEMBL alias seen in the source rows; otherwise a new Drug node is created (source = LAYER).

Relationship properties: layer = 'added curated drug-target links', sources, source_detail, action_labels,
push, mechanism, evidence, n_sources, drug_stage, release, added.

Dry run by default; --apply writes in one transaction and verifies counts; --undo removes exactly this layer.

    python code/write_new_links_to_graph.py
    python code/write_new_links_to_graph.py --apply
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, counts, driver, now, write_log, write_session  # noqa: E402

LAYER = "added curated drug-target links"
INP = ROOT / "results/graph_build/step0_4_new_links"
OUT = INP / "graph_write"


def text(v):
    return "" if pd.isna(v) else str(v)


def load_plan(sess):
    m = pd.read_csv(INP / "new_links_medicines_curated.csv")
    if m.empty or not m.evidence.eq("curated drug-target").all() or m.duplicated(["Drug_ID", "Gene_ID"]).any():
        raise SystemExit("expected unique curated drug-target pairs")
    raw = pd.read_csv(INP / "raw_rows_all_sources.csv.gz", low_memory=False, usecols=["Drug_key", "chembl"])
    aliases = raw.dropna().groupby("Drug_key").chembl.agg(lambda s: sorted(set(s))).to_dict()
    rel = json.loads((INP / "numbers.json").read_text(encoding="utf-8"))["release"]
    release = "; ".join(f"{k} {v}" for k, v in rel.items())
    cand = [{"key": d, "ids": [d] + [a for a in aliases.get(d, []) if a != d]} for d in sorted(m.Drug_ID.unique())]
    matches = sess.run("""UNWIND $rows AS row UNWIND range(0, size(row.ids)-1) AS i
        MATCH (d:Drug {node_id: row.ids[i]}) WITH row.key AS key, i, d.node_id AS id ORDER BY i
        RETURN key, collect(id)[0] AS id""", rows=cand).data()
    node_of = {r["key"]: r["id"] for r in matches}
    genes = set(sess.run("UNWIND $ids AS id MATCH (g:Gene {ensembl_id: id}) RETURN g.ensembl_id",
                         ids=sorted(m.Gene_ID.unique())).value())
    if genes != set(m.Gene_ID):
        raise SystemExit("some input genes are absent from the graph")
    new_nodes, rows = {}, []
    for r in m.itertuples():
        drug = node_of.get(r.Drug_ID, r.Drug_ID)
        if r.Drug_ID not in node_of:
            props = {"node_id": drug, "name": text(r.Drug), "drug_stage": text(r.drug_stage), "source": LAYER,
                     "release": release}
            if pd.notna(r.max_phase):
                props["max_phase"] = float(r.max_phase)
                props["max_phase_source"] = "step 0.4 sources (ChEMBL / Open Targets)"
            if pd.notna(r.molecule_type):
                props["molecule_type"] = str(r.molecule_type)
            new_nodes.setdefault(drug, props)          # first row supplies node metadata
        rows.append({"drug": drug, "gene": r.Gene_ID, "props": {
            "action_labels": text(r.action_labels), "push": text(r.push), "sources": text(r.sources),
            "source_detail": text(r.source_detail), "evidence": text(r.evidence), "n_sources": int(r.n_sources),
            "drug_stage": text(r.drug_stage), "mechanism": text(r.mechanism_text), "layer": LAYER,
            "release": release}})
    if len({(r["drug"], r["gene"]) for r in rows}) != len(rows):
        raise SystemExit("drug alias mapping collapses distinct pairs; resolve before loading")
    return {"nodes": list(new_nodes.values()), "links": rows, "release": release,
            "matched_existing_drugs": len(node_of)}


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--undo", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    drv = driver()
    with write_session(drv) as s:
        before = counts(s)
        if args.undo:
            def tx_undo(tx):
                n = tx.run("MATCH ()-[r:ACTS_ON {layer: $l}]->() DELETE r RETURN count(r) AS n", l=LAYER).single()["n"]
                k = tx.run("MATCH (d:Drug {source: $l}) WHERE NOT (d)--() DELETE d RETURN count(d) AS n",
                           l=LAYER).single()["n"]
                return n, k
            n, k = s.execute_write(tx_undo)
            write_log(OUT / "undo_log.json", {"undone": now(), "links_removed": n, "drug_nodes_removed": k,
                                              "before": before, "after": counts(s)})
            print(f"removed {n} links and {k} drug nodes")
            return
        if s.run("MATCH ()-[r:ACTS_ON {layer: $l}]->() RETURN count(r)", l=LAYER).single()[0]:
            raise SystemExit("this layer already exists; undo it first")
        plan = load_plan(s)
        print(f"planned: {len(plan['links']):,} links; {len(plan['nodes']):,} new Drug nodes; "
              f"{plan['matched_existing_drugs']:,} drug keys matched existing Drug nodes")
        pd.DataFrame([{"drug": r["drug"], "gene": r["gene"], **r["props"]} for r in plan["links"]]).to_csv(
            OUT / "links_planned.csv", index=False)
        pd.DataFrame(plan["nodes"]).to_csv(OUT / "drug_nodes_planned.csv", index=False)
        if not args.apply:
            print("dry run: nothing written")
            return
        added = now()
        (OUT / "undo_new_links.cypher").write_text(
            f"// Undo step 0.4 (or run: python code/write_new_links_to_graph.py --undo)\n"
            f"MATCH ()-[r:ACTS_ON {{layer: '{LAYER}'}}]->() DELETE r;\n"
            f"MATCH (d:Drug {{source: '{LAYER}'}}) WHERE NOT (d)--() DELETE d;\n", encoding="utf-8")

        def tx_write(tx):
            for i in range(0, len(plan["nodes"]), 1000):
                tx.run("UNWIND $rows AS row CREATE (d:Drug) SET d = row, d.added = $a",
                       rows=plan["nodes"][i:i + 1000], a=added).consume()
            n = 0
            for i in range(0, len(plan["links"]), 1000):
                n += tx.run("""UNWIND $rows AS row MATCH (d:Drug {node_id: row.drug})
                    MATCH (g:Gene {ensembl_id: row.gene}) CREATE (d)-[r:ACTS_ON]->(g)
                    SET r = row.props, r.added = $a RETURN count(r) AS n""",
                            rows=plan["links"][i:i + 1000], a=added).single()["n"]
            if n != len(plan["links"]):
                raise RuntimeError(f"wrote {n} of {len(plan['links'])} links; rolled back")
        s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    ok = (after["nodes"].get("Drug", 0) - before["nodes"].get("Drug", 0) == len(plan["nodes"])
          and after["relationships"].get("ACTS_ON", 0) - before["relationships"].get("ACTS_ON", 0) == len(plan["links"]))
    write_log(OUT / "write_log.json", {"step": "0.4 new curated drug-gene links", "added": added, "layer": LAYER,
                                       "links_written": len(plan["links"]), "drug_nodes_created": len(plan["nodes"]),
                                       "release": plan["release"], "before": before, "after": after,
                                       "count_check_passed": ok, "undo": "undo_new_links.cypher"})
    if not ok:
        raise SystemExit("count check failed")
    print(f"written: {len(plan['links']):,} links, {len(plan['nodes']):,} Drug nodes; count check passed")


if __name__ == "__main__":
    main()
