"""Step 0.6: add KEGG legacy and Hallmark pathways (MSigDB 2024.1.Hs GMT files) as Pathway nodes + IN_PATHWAY.

Parsing: node ids 'KEGG:<term>', readable names, members matched by gene
symbol; no expected counts are hard-coded and nothing is linked to NAFLD here (that is step 0.7).
Hallmark sets get node ids 'HALLMARK:<term>'.

Gene matching: a GMT symbol is matched to the Gene node(s) with that `symbol`. If a symbol belongs to more
than one Gene node, only protein-coding nodes are used; if still ambiguous, the symbol is skipped and listed.

Dry run by default; --apply writes in one transaction and checks counts. --undo removes the layer.

    python code/load_pathways.py
    python code/load_pathways.py --apply
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, counts, driver, now, write_log, write_session  # noqa: E402

OUT = ROOT / "results/graph_build/step0_6_pathways"
GMTS = [("KEGG_legacy", "KEGG:", "KEGG_", ROOT / "data/msigdb/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"),
        ("Hallmark", "HALLMARK:", "HALLMARK_", ROOT / "data/msigdb/h.all.v2024.1.Hs.symbols.gmt")]
ACRONYMS = {"ABC", "TGF", "ECM", "DNA", "RNA", "MAPK", "JAK", "STAT", "NOD", "RIG", "PPAR", "MTOR", "ERBB", "VEGF",
            "GNRH", "WNT", "NK", "T", "B", "II", "UV", "KRAS", "TNFA", "NFKB", "MYC", "E2F", "G2M", "UPR", "EMT",
            "SNARE", "GPI", "TCA", "FC", "HIF", "HCM", "ARVC", "SHH", "IL6", "IL2", "PI3K", "AKT", "P53"}


def readable(term, prefix):
    words = term[len(prefix):].split("_")
    out = [w if w in ACRONYMS or any(c.isdigit() for c in w) else w.lower() for w in words]
    s = " ".join(out)
    return s[0].upper() + s[1:]


def build_plan(sess):
    genes = pd.DataFrame(sess.run("MATCH (g:Gene) WHERE g.symbol IS NOT NULL "
                                  "RETURN g.node_id AS id, g.symbol AS symbol, g.biotype AS biotype").data())
    by_sym = {}
    for sym, grp in genes.groupby("symbol"):
        ids = list(grp.id)
        if len(ids) > 1:
            pc = list(grp[grp.biotype == "protein_coding"].id)
            ids = pc if len(pc) == 1 else ([] if len(pc) != 1 else pc)
        by_sym[sym] = ids
    nodes, members, unmatched, ambiguous = [], [], set(), set()
    for source, idp, prefix, path in GMTS:
        for line in path.read_text(encoding="utf-8").splitlines():
            parts = line.rstrip("\n").split("\t")
            term, gs = parts[0], sorted(set(parts[2:]) - {""})
            if len(parts) < 3 or not term.startswith(prefix):
                raise SystemExit(f"{path.name}: unexpected line {term}")
            node = idp + term
            nodes.append({"node_id": node, "name": readable(term, prefix), "msigdb_id": term, "source": source,
                          "msigdb_version": "2024.1.Hs", "n_genes": len(gs)})
            for g in gs:
                ids = by_sym.get(g)
                if ids is None:
                    unmatched.add(g)
                elif not ids:
                    ambiguous.add(g)
                for gid in ids or []:
                    members.append({"gene": gid, "pathway": node, "source": source})
    return nodes, members, sorted(unmatched), sorted(ambiguous)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--undo", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    drv = driver()
    with write_session(drv) as s:
        if args.undo:
            n = s.run("MATCH (p:Pathway) WHERE p.source IN ['KEGG_legacy','Hallmark'] DETACH DELETE p RETURN count(p)").single()[0]
            print("removed", n, "pathways")
            return
        if s.run("MATCH (p:Pathway) WHERE p.source IN ['KEGG_legacy','Hallmark'] RETURN count(p)").single()[0]:
            raise SystemExit("KEGG/Hallmark pathways already loaded; undo first")
        nodes, members, unmatched, ambiguous = build_plan(s)
        summ = pd.DataFrame(nodes).groupby("source").agg(pathways=("node_id", "size"), gmt_memberships=("n_genes", "sum"))
        summ["memberships_matched"] = pd.DataFrame(members).groupby("source").size()
        summ.to_csv(OUT / "pathway_summary.csv")
        pd.Series(unmatched, name="symbol").to_csv(OUT / "unmatched_symbols.csv", index=False)
        pd.Series(ambiguous, name="symbol").to_csv(OUT / "ambiguous_symbols.csv", index=False)
        print(summ.to_string())
        print(f"GMT symbols not found among Gene nodes: {len(unmatched)}; ambiguous (skipped): {len(ambiguous)}")
        if not args.apply:
            print("dry run: nothing written")
            return
        before = counts(s)
        added = now()
        (OUT / "undo_pathways.cypher").write_text(
            "// Undo step 0.6\nMATCH (p:Pathway) WHERE p.source IN ['KEGG_legacy','Hallmark'] DETACH DELETE p;\n",
            encoding="utf-8")

        def tx_write(tx):
            tx.run("UNWIND $rows AS row CREATE (p:Pathway) SET p = row, p.added = $a", rows=nodes, a=added).consume()
            n = 0
            for i in range(0, len(members), 5000):
                n += tx.run("""UNWIND $rows AS row MATCH (g:Gene {node_id: row.gene}) MATCH (p:Pathway {node_id: row.pathway})
                    CREATE (g)-[r:IN_PATHWAY]->(p) SET r.source = row.source RETURN count(r) AS n""",
                            rows=members[i:i + 5000]).single()["n"]
            if n != len(members):
                raise RuntimeError("membership count mismatch; rolled back")
        s.execute_write(tx_write)
        after = counts(s)
    drv.close()
    ok = (after["nodes"]["Pathway"] - before["nodes"]["Pathway"] == len(nodes)
          and after["relationships"]["IN_PATHWAY"] - before["relationships"]["IN_PATHWAY"] == len(members))
    write_log(OUT / "write_log.json", {"step": "0.6 KEGG legacy + Hallmark pathways", "added": added,
                                       "pathways": len(nodes), "memberships": len(members), "before": before,
                                       "after": after, "count_check_passed": ok, "undo": "undo_pathways.cypher"})
    if not ok:
        raise SystemExit("count check failed")
    print(f"written: {len(nodes)} pathways, {len(members)} memberships; count check passed")


if __name__ == "__main__":
    main()
