"""Step 0.2: fill drug metadata OptimusKG leaves empty, from ChEMBL (queried now).

Only drugs with a gap are
asked (no SMILES, no molecule type, no max phase, or no name, i.e. OptimusKG stored the accession as
the name), only ChEMBL accessions can be asked, and OptimusKG values always take precedence.

Written to Drug nodes only where the property is empty:
  smiles, molecule_type, max_phase          (+ <field>_source = 'ChEMBL <release>')
  recovered_name, recovered_name_source     (ChEMBL pref_name for drugs OptimusKG left unnamed)

Raw ChEMBL responses: data/raw/chembl_gaps/. Outputs: results/graph_build/step0_2_chembl_gaps/.

    python code/fetch_chembl_gaps.py           # fetch + dry run
    python code/fetch_chembl_gaps.py --apply
"""
import argparse
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, counts, driver, now, read_session, write_log, write_session  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data"
RAW = ROOT / "data/raw/chembl_gaps"
OUT = ROOT / "results/graph_build/step0_2_chembl_gaps"
BATCH = 20
FIELDS = ["smiles", "molecule_type", "max_phase"]


def get_json(url, cache):
    f = RAW / cache
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    last = None
    for attempt in range(5):
        try:
            js = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}),
                                                  timeout=120))
            f.write_text(json.dumps(js), encoding="utf-8")
            time.sleep(0.1)
            return js
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"ChEMBL request failed: {url[:120]} {last}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    drv = driver()
    with read_session(drv) as s:
        m = pd.DataFrame(s.run("""MATCH (d:Drug) RETURN d.node_id AS Drug_ID, d.name AS name,
            d.named_in_optimuskg AS named, d.smiles AS smiles, d.molecule_type AS molecule_type,
            d.max_phase AS max_phase""").data())
    need = m[m.smiles.isna() | m.molecule_type.isna() | m.max_phase.isna() | ~m.named.fillna(False).astype(bool)]
    ids = sorted(i for i in need.Drug_ID if str(i).startswith("CHEMBL"))
    print(f"Drug nodes {len(m):,}; with a gap {len(need):,}; ChEMBL accessions to ask {len(ids):,}")

    status = get_json(f"{API}/status.json", f"status_{now()[:10]}.json")
    release = f"{status.get('chembl_db_version')} ({status.get('chembl_release_date')})"
    rows = []

    def batch(i):
        return get_json(f"{API}/molecule.json?molecule_chembl_id__in={','.join(ids[i:i + BATCH])}&limit=100",
                        f"molecule_{i}.json")
    with ThreadPoolExecutor(6) as pool:
        pages = list(pool.map(batch, range(0, len(ids), BATCH)))
    for i, js in zip(range(0, len(ids), BATCH), pages):
        for x in js.get("molecules", []):
            st = x.get("molecule_structures") or {}
            rows.append({"Drug_ID": x.get("molecule_chembl_id"), "chembl_pref_name": x.get("pref_name"),
                         "chembl_molecule_type": x.get("molecule_type"), "chembl_max_phase": x.get("max_phase"),
                         "chembl_canonical_smiles": st.get("canonical_smiles")})
    c = pd.DataFrame(rows, columns=["Drug_ID", "chembl_pref_name", "chembl_molecule_type", "chembl_max_phase",
                                    "chembl_canonical_smiles"]).drop_duplicates("Drug_ID")
    c.to_csv(OUT / "chembl_gap_fill.csv", index=False)

    x = need.merge(c, on="Drug_ID", how="left")
    plan = []
    for r in x.itertuples():
        props = {}
        if pd.isna(r.smiles) and pd.notna(r.chembl_canonical_smiles):
            props["smiles"] = r.chembl_canonical_smiles
        if pd.isna(r.molecule_type) and pd.notna(r.chembl_molecule_type):
            props["molecule_type"] = r.chembl_molecule_type
        if pd.isna(r.max_phase) and pd.notna(r.chembl_max_phase):
            props["max_phase"] = float(r.chembl_max_phase)
        for f in list(props):
            props[f + "_source"] = f"ChEMBL {release}"
        if not (r.named is True) and pd.notna(r.chembl_pref_name):
            props["recovered_name"] = r.chembl_pref_name
            props["recovered_name_source"] = f"ChEMBL {release}"
        if props:
            plan.append({"id": r.Drug_ID, "props": props})
    filled = {f: sum(f in p["props"] for p in plan) for f in FIELDS + ["recovered_name"]}
    pd.DataFrame([{"Drug_ID": p["id"], **p["props"]} for p in plan]).to_csv(OUT / "planned_writes.csv", index=False)
    summary = {"chembl_release": release, "drug_nodes": len(m), "drugs_with_gap": len(need),
               "chembl_accessions_asked": len(ids), "chembl_records_returned": len(c),
               "gap_by_field_before": {"smiles": int(m.smiles.isna().sum()), "molecule_type": int(m.molecule_type.isna().sum()),
                                       "max_phase": int(m.max_phase.isna().sum()),
                                       "unnamed": int((~m.named.fillna(False).astype(bool)).sum())},
               "filled_by_chembl": filled, "drugs_written": len(plan)}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if not args.apply:
        print("dry run: nothing written")
        return

    allp = sorted({k for p in plan for k in p["props"]})
    (OUT / "undo_chembl_gaps.cypher").write_text(
        "// Undo step 0.2: remove only the properties ChEMBL filled (OptimusKG values were never overwritten)\n"
        + "".join(f"MATCH (d:Drug) WHERE d.{f}_source STARTS WITH 'ChEMBL' REMOVE d.{f}, d.{f}_source;\n"
                  for f in FIELDS + ["recovered_name"]), encoding="utf-8")
    with write_session(drv) as s:
        before = counts(s)
        started = now()

        def tx_write(tx):
            n = tx.run("""UNWIND $rows AS r MATCH (d:Drug {node_id: r.id})
                WITH d, r, [k IN keys(r.props) WHERE NOT k ENDS WITH '_source' AND k <> 'recovered_name'
                            AND d[k] IS NOT NULL] AS clash
                WHERE size(clash) = 0
                SET d += r.props RETURN count(d) AS n""", rows=plan).single()["n"]
            if n != len(plan):
                raise RuntimeError(f"wrote {n} of {len(plan)}; some OptimusKG value would be overwritten. Rolled back.")
            return n
        n = s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    if after["nodes"].get("Drug") != before["nodes"].get("Drug"):
        raise SystemExit("Drug node count changed; this step only sets properties")
    write_log(OUT / "write_log.json", {"step": "0.2 ChEMBL drug-metadata gaps", "started": started, "finished": now(),
                                       "drugs_written": n, "properties": allp, "summary": summary,
                                       "before": before, "after": after, "undo": "undo_chembl_gaps.cypher"})
    print(f"written: {n} Drug nodes")


if __name__ == "__main__":
    main()
