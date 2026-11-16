"""Step 0.5: drug action directions for OptimusKG ACTS_ON links on the seed genes that have none.

A link has no direction when its OptimusKG action gives none (TARGET / ENZYME / TRANSPORTER / CARRIER and other
non-directional labels; same action_direction() rule as Task 2). Each such link is looked up in three
public resources, queried now:
  1. ChEMBL drug mechanisms (REST API; drug ChEMBL ID incl. parent compound + target gene symbol)
  2. IUPHAR/BPS Guide to Pharmacology (interactions.csv downloaded in step 0.4; ligand ChEMBL ID + Ensembl ID)
  3. DGIdb 5 (GraphQL API; gene symbol + drug ChEMBL ID or name)
A link gets a direction only if every source that gives one agrees ('down' or 'up'); otherwise
'conflicting' or 'still unknown'. Resolved links get new properties (never changing `action`):
  external_action_direction, external_action_label, external_action_sources, external_action_release,
  external_action_added

Raw downloads: data/raw/external_directions/. Outputs: results/graph_build/step0_5_external_directions/.

    python code/external_action_check.py            # lookup + dry run
    python code/external_action_check.py --apply
    python code/external_action_check.py --undo
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, counts, driver, now, read_session, write_log, write_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
OUT = ROOT / "results/graph_build/step0_5_external_directions"
RAW = ROOT / "data/raw/external_directions"
GTOP = ROOT / "data/raw/new_links"
PROPS = ["external_action_direction", "external_action_label", "external_action_sources",
         "external_action_release", "external_action_added"]

DOWN_WORDS = ("INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE", "INVERSE AGONIST", "INVERSE_AGONIST", "DEGRADER",
              "INHIBITION", "INHIBITORY", "ANTISENSE", "RNAI", "DISRUPTING", "SUPPRESSOR", "CHANNEL BLOCKER",
              "GATING INHIBITOR", "ALLOSTERIC ANTAGONIST")
UP_WORDS = ("AGONIST", "ACTIVATOR", "POSITIVE", "OPENER", "STABILISER", "STABILIZER", "INDUCER", "ACTIVATION",
            "ACTIVATING", "POTENTIATOR", "ENHANCER")
# Task 2 rule for an OptimusKG action
DOWN = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE_MODULATOR", "NEGATIVE_ALLOSTERIC_MODULATOR",
        "ALLOSTERIC_ANTAGONIST", "INVERSE_AGONIST", "DEGRADER"}
UP = {"AGONIST", "PARTIAL_AGONIST", "ACTIVATOR", "POSITIVE_MODULATOR", "POSITIVE_ALLOSTERIC_MODULATOR",
      "OPENER", "STABILISER"}
MIXED = {"MODULATOR", "BINDING_AGENT", "SUBSTRATE", "RELEASING_AGENT"}


def direction(label):
    s = str(label).upper().replace("_", " ")
    if not s or s in ("NAN", "NONE"):
        return None
    if any(w in s for w in DOWN_WORDS):
        return "down"
    if any(w in s for w in UP_WORDS):
        return "up"
    return None


def action_direction(a):
    a = str(a).upper()
    return "down" if a in DOWN else "up" if a in UP else "mixed" if a in MIXED else "unknown"


def edge_direction(action, assertions):
    d = action_direction(action)
    if isinstance(assertions, str) and assertions:
        if {"down", "up"} <= {action_direction(x.split(":", 1)[-1]) for x in assertions.split(";")}:
            return "mixed"
    return d


def get_json(url, cache_name, data=None):
    f = RAW / cache_name
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json", "Content-Type": "application/json"})
    last = None
    for attempt in range(5):
        try:
            js = json.load(urllib.request.urlopen(req, timeout=120))
            f.write_text(json.dumps(js), encoding="utf-8")
            time.sleep(0.2)
            return js
        except Exception as ex:  # noqa: BLE001
            last = ex
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"failed: {url[:120]} {last}")


def pick(df, keys, label):
    g = df.dropna(subset=["dir"]).groupby(keys).agg(dirs=("dir", lambda s: ";".join(sorted(set(s)))),
                                                  actions=("action", lambda s: "; ".join(sorted(set(map(str, s))))[:200]))
    g[f"{label}_dir"] = g.dirs.map(lambda s: "mixed" if ";" in s else s)
    return g[[f"{label}_dir"]].join(g.actions.rename(f"{label}_action"))


def lookup():
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    drv = driver()
    with read_session(drv) as s:
        e = pd.DataFrame(s.run("""MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids AND r.layer IS NULL
            RETURN d.node_id AS Drug_ID, coalesce(d.recovered_name, d.name) AS Drug, g.ensembl_id AS Gene_ID,
                   g.symbol AS Gene, r.action AS Action_raw, r.relation_assertions AS assertions,
                   r.external_action_direction AS already""", ids=list(seeds.ensembl_id)).data())
    drv.close()
    if e.already.notna().any():
        raise SystemExit("external directions already written; undo first")
    e["Action_direction"] = [edge_direction(a, s) for a, s in zip(e.Action_raw, e.assertions)]
    unk = e[e.Action_direction == "unknown"].copy()
    drugs = sorted(d for d in e.Drug_ID.unique() if str(d).startswith("CHEMBL"))
    print(f"OptimusKG links on seed genes {len(e):,}; without direction {len(unk):,} "
          f"({unk.Action_raw.value_counts().to_dict()}); drugs {len(drugs):,}")

    # 1. ChEMBL mechanisms
    status = get_json("https://www.ebi.ac.uk/chembl/api/data/status.json", "chembl_status.json")
    chembl_release = f"ChEMBL {status.get('chembl_db_version')}"
    parent = {}
    for i in range(0, len(drugs), 40):
        js = get_json("https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit=1000&only=molecule_chembl_id,molecule_hierarchy"
                      "&molecule_chembl_id__in=" + ",".join(drugs[i:i + 40]), f"chembl_molecule_{i}.json")
        for m in js["molecules"]:
            h = m.get("molecule_hierarchy") or {}
            parent[m["molecule_chembl_id"]] = h.get("parent_chembl_id") or m["molecule_chembl_id"]
    qids = sorted(set(drugs) | set(parent.values()))
    mech = []
    for i in range(0, len(qids), 40):
        b, offset = qids[i:i + 40], 0
        while True:
            js = get_json("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?limit=1000&offset=%d&molecule_chembl_id__in=%s"
                          % (offset, ",".join(b)), f"chembl_mech_{i}_{offset}.json")
            mech += js["mechanisms"]
            if not js["page_meta"]["next"]:
                break
            offset += 1000
    mech = pd.DataFrame(mech)
    targets = sorted(mech.target_chembl_id.dropna().unique()) if len(mech) else []
    tsym = {}
    for i in range(0, len(targets), 40):
        js = get_json("https://www.ebi.ac.uk/chembl/api/data/target.json?limit=1000&only=target_chembl_id,target_components"
                      "&target_chembl_id__in=" + ",".join(targets[i:i + 40]), f"chembl_target_{i}.json")
        for t in js["targets"]:
            tsym[t["target_chembl_id"]] = {sy["component_synonym"].upper() for c in t.get("target_components") or []
                                           for sy in c.get("target_component_synonyms") or []
                                           if sy.get("syn_type") == "GENE_SYMBOL"}
    child_of = {}
    for d, p in parent.items():
        child_of.setdefault(p, set()).add(d)
    rows = [(d, sym, m.action_type, m.mechanism_of_action) for m in mech.itertuples()
            for sym in tsym.get(m.target_chembl_id, ()) for d in {m.molecule_chembl_id} | child_of.get(m.molecule_chembl_id, set())]
    chembl = pd.DataFrame(rows, columns=["Drug_ID", "symbol", "action", "mechanism"]).drop_duplicates()
    chembl["dir"] = chembl.action.map(direction)
    print(f"ChEMBL: {len(mech):,} mechanisms, {len(chembl):,} drug-gene rows")

    # 2. Guide to Pharmacology (files downloaded in step 0.4)
    gfile = GTOP / "gtopdb_interactions.csv"
    gtop_version = gfile.read_text(encoding="utf-8").splitlines()[0].strip('"')
    gi = pd.read_csv(gfile, skiprows=1, low_memory=False)
    gm = pd.read_csv(GTOP / "gtopdb_ligand_id_mapping.csv", skiprows=1, low_memory=False)
    lig_col = next(c for c in gm.columns if c.lower() == "ligand id")
    chembl_col = next(c for c in gm.columns if c.lower() == "chembl id")
    gi = gi[gi["Target Species"].fillna("").str.contains("Human", case=False)].copy()
    gi["Drug_ID"] = gi["Ligand ID"].astype(str).map(dict(zip(gm[lig_col].astype(str), gm[chembl_col])))
    gi = gi.dropna(subset=["Drug_ID", "Target Ensembl Gene ID"])
    gi = gi.assign(Gene_ID=gi["Target Ensembl Gene ID"].astype(str).str.split("|")).explode("Gene_ID")
    gtop = pd.DataFrame({"Drug_ID": gi.Drug_ID, "Gene_ID": gi.Gene_ID.str.strip(),
                         "action": gi["Type"].fillna("") + " / " + gi["Action"].fillna("")})
    gtop["dir"] = [direction(t) or direction(a) for t, a in zip(gi["Type"], gi["Action"])]
    print(f"GtoPdb ({gtop_version}): {len(gtop):,} human ligand-target rows")

    # 3. DGIdb GraphQL
    genes_u = sorted(set(e.Gene.dropna()))
    drows = []
    for i in range(0, len(genes_u), 25):
        q = ('{ genes(names: %s) { nodes { name interactions { drug { name conceptId } interactionTypes { type directionality } '
             'sources { sourceDbName } } } } }') % json.dumps(genes_u[i:i + 25])
        js = get_json("https://dgidb.org/api/graphql", f"dgidb_{i}.json", data=json.dumps({"query": q}).encode())
        for gnode in js["data"]["genes"]["nodes"]:
            for it in gnode["interactions"]:
                types = it["interactionTypes"] or []
                dirs = {("down" if t["directionality"] == "INHIBITORY" else "up" if t["directionality"] == "ACTIVATING" else None)
                        for t in types} - {None}
                drows.append((gnode["name"].upper(), it["drug"]["name"].upper(), it["drug"]["conceptId"],
                              ";".join(sorted(t["type"] for t in types)), ";".join(sorted(dirs)),
                              ";".join(sorted(s_["sourceDbName"] for s_ in it["sources"]))))
    dg = pd.DataFrame(drows, columns=["symbol", "drug_name", "concept", "action", "dirs", "sources"])
    dg["dir"] = dg.dirs.map(lambda s: None if not s else ("mixed" if ";" in s else s))
    dg["Drug_ID"] = dg.concept.str.extract(r"chembl:(CHEMBL\d+)", expand=False)
    print(f"DGIdb: {len(dg):,} gene-drug interactions for {len(genes_u):,} genes")

    def attach(df):
        df = df.copy()
        df["symbol"] = df.Gene.str.upper()
        df["drug_upper"] = df.Drug.str.upper()
        df = df.join(pick(chembl, ["Drug_ID", "symbol"], "chembl"), on=["Drug_ID", "symbol"])
        df = df.join(pick(gtop, ["Drug_ID", "Gene_ID"], "gtopdb"), on=["Drug_ID", "Gene_ID"])
        df = df.join(pick(dg.dropna(subset=["Drug_ID"]), ["Drug_ID", "symbol"], "dgidb_id"), on=["Drug_ID", "symbol"])
        df = df.join(pick(dg, ["drug_name", "symbol"], "dgidb_nm").rename_axis(["drug_upper", "symbol"]),
                     on=["drug_upper", "symbol"])
        df["dgidb_dir"] = df.dgidb_id_dir.fillna(df.dgidb_nm_dir)
        df["dgidb_action"] = df.dgidb_id_action.fillna(df.dgidb_nm_action)
        return df

    u = attach(unk)

    def combine(r):
        ds = {d for d in (r.chembl_dir, r.gtopdb_dir, r.dgidb_dir) if isinstance(d, str)}
        if not ds:
            return "still unknown"
        return "down" if ds == {"down"} else "up" if ds == {"up"} else "conflicting"

    u["external_direction"] = u.apply(combine, axis=1) if len(u) else []
    u["n_sources_agreeing"] = u.apply(lambda r: sum(isinstance(x, str) and x == r.external_direction
                                                    for x in (r.chembl_dir, r.gtopdb_dir, r.dgidb_dir)), axis=1) if len(u) else []
    u["graph_action_type"] = np.where(u.Action_raw == "TARGET", "TARGET (binding, no direction)",
                                      np.where(u.Action_raw.isin(["ENZYME", "TRANSPORTER", "CARRIER"]),
                                               "ENZYME/TRANSPORTER/CARRIER", "other non-directional label"))
    u = u.merge(seeds[["ensembl_id", "main_stage", "direction"]].rename(
        columns={"ensembl_id": "Gene_ID", "direction": "gene_direction"}), on="Gene_ID", how="left")
    keep = ["Drug_ID", "Drug", "Gene_ID", "Gene", "Action_raw", "graph_action_type", "main_stage", "gene_direction",
            "chembl_dir", "chembl_action", "gtopdb_dir", "gtopdb_action", "dgidb_dir", "dgidb_action",
            "external_direction", "n_sources_agreeing"]
    u[keep].to_csv(OUT / "unknown_links_external_lookup.csv", index=False)

    k = attach(e[e.Action_direction.isin(["up", "down"])][["Drug_ID", "Drug", "Gene_ID", "Gene", "Action_direction"]])
    ctrl = []
    for lab, col in [("ChEMBL mechanisms", "chembl_dir"), ("Guide to Pharmacology", "gtopdb_dir"), ("DGIdb", "dgidb_dir")]:
        f_ = k[col].isin(["up", "down"])
        a_ = f_ & (k[col] == k.Action_direction)
        ctrl.append((lab, len(k), int(f_.sum()), int(a_.sum()), round(100 * a_.sum() / max(f_.sum(), 1), 1)))
    ctrl = pd.DataFrame(ctrl, columns=["source", "known-direction links", "found", "agree with graph", "% agree"])
    ctrl.to_csv(OUT / "positive_control.csv", index=False)

    release = (f"{chembl_release} mechanisms (REST API); {gtop_version}; DGIdb 5 GraphQL API; queried {now()[:10]}")
    resolved = u[u.external_direction.isin(["down", "up"])]
    rows = [{"drug": d, "gene": g, "direction": p,
             "label": " | ".join(str(x) for x in (ca, ga, da) if isinstance(x, str))[:300],
             "sources": ";".join(n for n, x in (("ChEMBL", cd), ("GtoPdb", gd), ("DGIdb", dd)) if x == p)}
            for d, g, p, ca, ga, da, cd, gd, dd in zip(resolved.Drug_ID, resolved.Gene_ID, resolved.external_direction,
                                                       resolved.chembl_action, resolved.gtopdb_action, resolved.dgidb_action,
                                                       resolved.chembl_dir, resolved.gtopdb_dir, resolved.dgidb_dir)]
    summary = {"run_utc": now(), "release": release, "optimuskg_links_on_seed_genes": len(e),
               "known_direction_before": int(e.Action_direction.isin(["up", "down"]).sum()),
               "mixed_before": int((e.Action_direction == "mixed").sum()),
               "no_direction": len(u), "no_direction_by_action": u.Action_raw.value_counts().to_dict(),
               "result": u.external_direction.value_counts().to_dict(),
               "by_graph_label": pd.crosstab(u.graph_action_type, u.external_direction).to_dict(),
               "links_to_write": len(rows), "positive_control": ctrl.to_dict(orient="records")}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))
    return rows, release


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--undo", action="store_true")
    args = ap.parse_args()
    if args.undo:
        drv = driver()
        with write_session(drv) as s:
            n = s.run("MATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL REMOVE "
                      + ", ".join(f"r.{p}" for p in PROPS) + " RETURN count(r)").single()[0]
        print("removed external properties from", n, "relationships")
        return
    rows, release = lookup()
    if not args.apply:
        print("dry run: nothing written")
        return
    drv = driver()
    with write_session(drv) as s:
        before = counts(s)
        bk = s.run("""UNWIND $rows AS row MATCH (:Drug {node_id: row.drug})-[r:ACTS_ON]->(:Gene {ensembl_id: row.gene})
            WHERE r.layer IS NULL RETURN row.drug AS drug, row.gene AS gene, properties(r) AS p""", rows=rows).data()
        if len(bk) != len(rows):
            raise SystemExit(f"{len(bk)} relationships matched for {len(rows)} pairs; expected exactly one each")
        pd.DataFrame([{"drug": b["drug"], "gene": b["gene"], "properties": json.dumps(b["p"], default=str)} for b in bk]).to_csv(
            OUT / "backup_ACTS_ON_properties_before_write.csv", index=False)
        (OUT / "undo_external_properties.cypher").write_text(
            "// Removes ONLY the properties added in step 0.5\nMATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL\n"
            "REMOVE " + ", ".join(f"r.{p}" for p in PROPS) + "\nRETURN count(r) AS removed;\n", encoding="utf-8")
        added = now()

        def tx_write(tx):
            n = tx.run("""UNWIND $rows AS row MATCH (:Drug {node_id: row.drug})-[r:ACTS_ON]->(:Gene {ensembl_id: row.gene})
                WHERE r.layer IS NULL
                SET r.external_action_direction = row.direction, r.external_action_label = row.label,
                    r.external_action_sources = row.sources, r.external_action_release = $rel,
                    r.external_action_added = $a RETURN count(r) AS n""", rows=rows, rel=release, a=added).single()["n"]
            if n != len(rows):
                raise RuntimeError("count mismatch; rolled back")
            return n
        n = s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    write_log(OUT / "write_log.json", {"step": "0.5 external directions", "added": added, "written": n,
                                       "release": release, "before": before, "after": after,
                                       "backup": "backup_ACTS_ON_properties_before_write.csv",
                                       "undo": "undo_external_properties.cypher"})
    print("written:", n)


if __name__ == "__main__":
    main()
