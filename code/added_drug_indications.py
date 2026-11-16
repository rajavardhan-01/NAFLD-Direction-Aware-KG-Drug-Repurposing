"""Step 0.8: ChEMBL indications for the Drug nodes added in step 0.4, written as TREATS.

Method:
  * added drugs = Drug nodes with source = 'added curated drug-target links' (they have no disease links)
  * ChEMBL ID = the node id if it is a ChEMBL accession, else an exact pref-name / synonym match
  * indications fetched for the molecule and as parent molecule (salt forms), queried now
  * one TREATS per drug-disease pair (highest max_phase_for_ind), only for diseases already in the graph
    (EFO id with ':' -> '_'), skipping pairs that already have TREATS; no node is created

Outputs: results/graph_build/step0_8_added_drug_indications/. Dry run by default; --apply; --undo.

    python code/added_drug_indications.py
    python code/added_drug_indications.py --apply
"""
import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, counts, driver, now, read_session, write_log, write_session  # noqa: E402

OUT = ROOT / "results/graph_build/step0_8_added_drug_indications"
ADDED_LAYER = "added curated drug-target links"
LAYER = "ChEMBL indication (added drugs)"
API = "https://www.ebi.ac.uk/chembl/api/data"
BATCH = 40
FIELDS = ["molecule_chembl_id", "parent_molecule_chembl_id", "efo_id", "efo_term", "mesh_id", "mesh_heading",
          "max_phase_for_ind", "drugind_id"]


def get(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=120) as r:
                return json.load(r)
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def indications(ids, field):
    rows, url = [], f"{API}/drug_indication.json?" + urllib.parse.urlencode({f"{field}__in": ",".join(ids), "limit": 1000})
    while url:
        d = get(url)
        rows += [{k: r.get(k) for k in FIELDS} for r in d["drug_indications"]]
        nxt = d["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def chembl_by_name(name):
    for field, how in [("pref_name__iexact", "preferred name"), ("molecule_synonyms__molecule_synonym__iexact", "synonym")]:
        d = get(f"{API}/molecule.json?" + urllib.parse.urlencode({field: name, "limit": 5}))
        hits = {m["molecule_chembl_id"] for m in d["molecules"]}
        if len(hits) == 1:
            return hits.pop(), how
        if len(hits) > 1:
            return None, f"ambiguous ({how}: {';'.join(sorted(hits))})"
    return None, "not found"


def fetch(drv):
    with read_session(drv) as s:
        added = pd.DataFrame(s.run("MATCH (d:Drug {source: $l}) RETURN d.node_id AS Drug_ID, d.name AS Drug, "
                                   "d.drug_stage AS drug_stage", l=ADDED_LAYER).data())
    print(f"added Drug nodes: {len(added)}")
    if added.empty:
        raise SystemExit("no added drugs; run step 0.4 first")
    is_chembl = added.Drug_ID.str.upper().str.startswith("CHEMBL")
    added["chembl_id"] = added.Drug_ID.where(is_chembl).str.upper()
    added["chembl_match"] = is_chembl.map({True: "graph ID", False: None})
    todo = list(added.index[~is_chembl])
    with ThreadPoolExecutor(max_workers=8) as pool:
        for i, (cid, how) in zip(todo, pool.map(lambda j: chembl_by_name(added.at[j, "Drug"]), todo)):
            added.at[i, "chembl_id"], added.at[i, "chembl_match"] = cid, how
    added.to_csv(OUT / "added_drug_chembl_ids.csv", index=False)
    lookup = added.dropna(subset=["chembl_id"])
    ids = sorted(set(lookup.chembl_id))
    status = get(f"{API}/status.json")
    rows = []
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        rows += indications(chunk, "molecule_chembl_id") + indications(chunk, "parent_molecule_chembl_id")
    ind = pd.DataFrame(rows, columns=FIELDS).drop_duplicates("drugind_id")
    ind["chembl_id"] = ind.molecule_chembl_id.where(ind.molecule_chembl_id.isin(ids), ind.parent_molecule_chembl_id)
    ind = ind[ind.chembl_id.isin(ids)].merge(lookup[["chembl_id", "Drug_ID", "Drug", "chembl_match"]], on="chembl_id")
    ind["max_phase_for_ind"] = pd.to_numeric(ind.max_phase_for_ind, errors="coerce")
    ind["efo_id_graph"] = ind.efo_id.str.replace(":", "_", regex=False)
    ind["chembl_release"] = status.get("chembl_db_version")
    ind["chembl_release_date"] = status.get("chembl_release_date")
    ind["retrieved"] = now()[:10]
    ind.sort_values(["Drug", "max_phase_for_ind"], ascending=[True, False]).to_csv(
        OUT / "added_drug_indications_chembl.csv", index=False)
    meta = {"chembl_release": status.get("chembl_db_version"), "chembl_release_date": status.get("chembl_release_date"),
            "retrieved": now()[:10], "added_drugs": len(added),
            "chembl_match": added.chembl_match.fillna("none").value_counts().to_dict(),
            "added_drugs_with_chembl_id": int(added.chembl_id.notna().sum()),
            "added_drugs_with_indication": int(ind.Drug_ID.nunique()), "indication_rows": len(ind)}
    (OUT / "added_drug_indications_chembl.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))
    return ind, meta


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--undo", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    drv = driver()
    if args.undo:
        with write_session(drv) as s:
            n = s.run("MATCH ()-[r:TREATS {layer: $l}]->() DELETE r RETURN count(r)", l=LAYER).single()[0]
        print("removed", n)
        return
    ind, meta = fetch(drv)
    pairs = (ind.groupby(["Drug_ID", "efo_id_graph"])
             .agg(Drug=("Drug", "first"), max_phase_for_ind=("max_phase_for_ind", "max"),
                  efo_terms=("efo_term", lambda s: "; ".join(sorted(set(s.dropna())))),
                  mesh_headings=("mesh_heading", lambda s: "; ".join(sorted(set(s.dropna())))),
                  drugind_ids=("drugind_id", lambda s: ";".join(str(int(x)) for x in sorted(set(s.dropna())))))
             .reset_index())
    with read_session(drv) as s:
        found = pd.DataFrame(s.run("""UNWIND $rows AS r
            OPTIONAL MATCH (d:Drug {node_id: r.drug}) OPTIONAL MATCH (x:Disease {node_id: r.dis})
            OPTIONAL MATCH (d)-[t:TREATS]->(x)
            RETURN r.drug AS Drug_ID, r.dis AS efo_id_graph, d IS NOT NULL AS drug_in_graph,
                   x IS NOT NULL AS disease_in_graph, count(t) AS existing_treats""",
                                   rows=[{"drug": a, "dis": b} for a, b in zip(pairs.Drug_ID, pairs.efo_id_graph)]).data())
        if s.run("MATCH ()-[r:TREATS {layer: $l}]->() RETURN count(r)", l=LAYER).single()[0]:
            raise SystemExit("this layer already exists; undo first")
    pairs = pairs.merge(found, on=["Drug_ID", "efo_id_graph"], how="left")
    pairs["status"] = "write"
    pairs.loc[~pairs.disease_in_graph.astype(bool), "status"] = "skip: disease not in graph"
    pairs.loc[~pairs.drug_in_graph.astype(bool), "status"] = "skip: drug not in graph"
    pairs.loc[pairs.existing_treats.astype(int) > 0, "status"] = "skip: TREATS already present"
    pairs.to_csv(OUT / "indication_pairs_to_write.csv", index=False)
    todo = pairs[pairs.status == "write"]
    print(pairs.status.value_counts().to_string())
    print(f"links to write: {len(todo):,} (drugs {todo.Drug_ID.nunique():,}, diseases {todo.efo_id_graph.nunique():,})")
    if not args.apply:
        print("dry run: nothing written")
        return
    added = now()
    rows = [{"drug": r.Drug_ID, "dis": r.efo_id_graph,
             "p": {"source": "ChEMBL indication", "layer": LAYER,
                   "max_phase_for_ind": None if pd.isna(r.max_phase_for_ind) else float(r.max_phase_for_ind),
                   "efo_terms": r.efo_terms, "mesh_headings": r.mesh_headings, "drugind_ids": r.drugind_ids,
                   "chembl_release": meta["chembl_release"], "chembl_release_date": meta["chembl_release_date"],
                   "retrieved": meta["retrieved"], "added": added}} for r in todo.itertuples()]
    with write_session(drv) as s:
        before = counts(s)
        (OUT / "undo_indications.cypher").write_text(f"MATCH ()-[r:TREATS {{layer: '{LAYER}'}}]->() DELETE r;\n",
                                                     encoding="utf-8")

        def tx_write(tx):
            n = tx.run("""UNWIND $rows AS r MATCH (d:Drug {node_id: r.drug}) MATCH (x:Disease {node_id: r.dis})
                CREATE (d)-[t:TREATS]->(x) SET t = r.p RETURN count(t) AS n""", rows=rows).single()["n"]
            if n != len(rows):
                raise RuntimeError("count mismatch; rolled back")
            return n
        n = s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    write_log(OUT / "write_log.json", {"step": "0.8 added-drug indications", "added": added, "links_written": n,
                                       "drugs": int(todo.Drug_ID.nunique()), "diseases": int(todo.efo_id_graph.nunique()),
                                       "skipped": pairs.status.value_counts().drop("write", errors="ignore").to_dict(),
                                       "meta": meta, "before": before, "after": after, "undo": "undo_indications.cypher"})
    print(f"written: {n} TREATS links")


if __name__ == "__main__":
    main()
