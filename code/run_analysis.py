"""Run the analyses (Tasks 1-4, direction filter, metapaths, ranking) against the paper2public graph.

Before running, checks that the graph has exactly this run's expected node and relationship counts
(results/graph_build/expected_graph_counts.json, written at the end of graph preparation). If the graph differs,
nothing runs. Read-only: no analysis step writes to Neo4j.

    python code/run_analysis.py --stage core        # Tasks 1-4 + direction filter
    python code/run_analysis.py --stage metapaths   # metapaths + ranking
    python code/run_analysis.py --stage all
    python code/run_analysis.py --write-expected    # record the current graph counts as expected (after graph build)
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, TARGET_DB, counts, driver, now, read_session  # noqa: E402

EXPECTED = ROOT / "results/graph_build/expected_graph_counts.json"
CORE = ["task1_seed_gene_table.py", "task2_drug_gene_edges.py", "task3_druggability.py", "task4_drug_coverage.py",
        "direction_consistent_drugs.py"]
METAPATHS = ["metapaths.py", "ranking.py"]


def graph_counts():
    drv = driver()
    with read_session(drv) as s:
        c = counts(s)
    drv.close()
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["core", "metapaths", "all"], default="all")
    ap.add_argument("--write-expected", action="store_true")
    args = ap.parse_args()
    c = graph_counts()
    if args.write_expected:
        EXPECTED.write_text(json.dumps({"database": TARGET_DB, "recorded": now(), **c}, indent=2), encoding="utf-8")
        print(f"recorded expected counts: {c['nodes_total']:,} nodes, {c['relationships_total']:,} relationships")
        return
    exp = json.loads(EXPECTED.read_text(encoding="utf-8"))
    diff = {k: (exp[k], c[k]) for k in ["nodes_total", "relationships_total", "nodes", "relationships"] if exp[k] != c[k]}
    if diff:
        raise SystemExit(f"graph does not match this run's expected counts; nothing run: {diff}")
    print(f"graph check passed: {c['nodes_total']:,} nodes, {c['relationships_total']:,} relationships")
    steps = (CORE if args.stage in ("core", "all") else []) + (METAPATHS if args.stage in ("metapaths", "all") else [])
    for sc in steps:
        t = time.time()
        r = subprocess.run([sys.executable, str(ROOT / "code" / sc)], cwd=ROOT)
        if r.returncode:
            raise SystemExit(f"{sc} failed (exit {r.returncode})")
        print(f"----- {sc} done in {time.time() - t:.0f} s", flush=True)


if __name__ == "__main__":
    main()
