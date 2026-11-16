"""Final audit (checks 1-4). Check 5 (rerun from scratch) is done by --rerun.

1. search code/ and README.md for the excluded inputs, paths and concepts (must be empty)
2. query paper2public: no expression nodes, no TRACKS_* relationships; list labels and relationship types
3. every Gene* in the route files is in the seed file (RunD: seed file or GWAS add-on file)
4. funnel numbers add up; direction groups sum to the candidate list; README equals a fresh build from the files
5. (--rerun) move the metapath/ranking outputs aside, rerun metapaths.py + ranking.py from scratch,
   and confirm the final ranked lists (and route files) are identical

Writes results/audit/audit_report.md and audit.json.

    python code/final_audit.py
    python code/final_audit.py --rerun
"""
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, TARGET_DB, counts, driver, now, read_session  # noqa: E402

OUT = ROOT / "results/audit"
RUNS = ROOT / "results/Metapath_runs"
RK = ROOT / "results/Ranking"
LOGS = ROOT / "logs"
FORBIDDEN = [r"liver_graph", r"graph-export-2026-11-05", r"rebuild_graph", r"prepare_graph_release", r"Paper1_HighConfidence",
             r"Paper1_Results", r"progression_genes", r"11_final_all_3845_genes", r"11_FINAL_shortlist", r"10_loso_all_eligible_runs",
             r"task2_gene_clusters", r"task4_all_significant", r"New_drug_gene_links", r"Task7_external_layer",
             r"indication_pairs_plan", r"graph_gap_properties", r"Tier [123]", r"[Ss]hortlist", r"23-gene", r"\bC[1-5]\b",
             r"\bF[0-4]\b", r"[Ll]argest_transition", r"6,222", r"3,845", r"14,794", r"Paper ?1", r"paper1", r"companion",
             r"old_reference", r"TRACKS_FIBROSIS", r"TRACKS_INFLAMMATION", r"Paper2_Results"]
ROUTE_FILES = ["{r}_BB.csv", "{r}_all_routes.csv"]
RANK_FILES = ["{r}_FINAL_ranked_drugs.csv", "{r}_strict_liver_treating_drugs.csv"]
RUN_NAMES = ["RunA_all_seed_genes", "RunB_histology_genes", "RunC_fibrosis_genes", "RunD_seed_plus_GWAS"]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def table_hash(p, sort_cols=None):
    d = pd.read_csv(p, low_memory=False)
    d = d.drop(columns=[c for c in d.columns if c in ("generated",)], errors="ignore")
    d = d.sort_values(sort_cols or list(d.columns)).reset_index(drop=True)
    return hashlib.sha256(d.to_csv(index=False).encode()).hexdigest()


def check1():
    hits = []
    files = sorted((ROOT / "code").glob("*.py")) + [ROOT / "README.md", ROOT / ".gitignore", ROOT / "requirements.txt"]
    pat = re.compile("|".join(FORBIDDEN))
    for f in files:
        if f.name == "final_audit.py":
            continue                       # this file holds the search terms themselves
        for i, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if pat.search(line):
                hits.append(f"{f.relative_to(ROOT)}:{i}: {line.strip()[:120]}")
    return {"files_searched": len(files) - 0, "patterns": FORBIDDEN, "hits": hits, "pass": not hits}


def check2():
    drv = driver()
    with read_session(drv) as s:
        c = counts(s)
        prop_hits = s.run("""MATCH (g:Gene) WHERE any(k IN keys(g) WHERE toLower(k) CONTAINS 'expression'
                             OR toLower(k) CONTAINS 'rho' OR toLower(k) CONTAINS 'fibrosis') RETURN count(g)""").single()[0]
    drv.close()
    tracks = [t for t in c["relationships"] if t.upper().startswith("TRACKS")]
    expr = [lab for lab in c["nodes"] if re.search("express|sample|study|tissue|anatomy", lab, re.I)]
    return {"database": TARGET_DB, "labels": c["nodes"], "relationship_types": c["relationships"],
            "tracks_relationship_types": tracks, "expression_like_labels": expr,
            "gene_nodes_with_expression_like_properties": prop_hits,
            "pass": not tracks and not expr and prop_hits == 0}


def check3():
    seed = set(pd.read_csv(ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv").ensembl_id)
    gwas = set(pd.read_csv(ROOT / "seed_genes/seed_genes_optional_GWAS_Catalog.csv").ensembl_id)
    out = {}
    for r in RUN_NAMES:
        allowed = seed | gwas if r.startswith("RunD") else seed
        g = set(pd.read_csv(RUNS / f"{r}_all_routes.csv").gene_star_id) | set(pd.read_csv(RUNS / f"{r}_BB.csv").gene_star_id)
        out[r] = {"gene_star_in_routes": len(g), "not_in_allowed_set": sorted(g - allowed)}
    return {"runs": out, "pass": all(not v["not_in_allowed_set"] for v in out.values())}


def check4():
    f = pd.read_csv(ROOT / "results/Task4/candidate_funnel.csv")
    steps_ok = all(f.drugs.iloc[i - 1] - f.drugs.iloc[i] == f.removed_at_step.iloc[i] for i in range(1, len(f)))
    mono = f.drugs.is_monotonic_decreasing
    main_ = pd.read_csv(ROOT / "results/Task4/repurposing_candidates_mechanism_annotated.csv")
    d = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_drugs_direction_summary.csv")
    only = pd.read_csv(ROOT / "results/Direction_consistent_drugs/only_consistent_drugs.csv")
    dir_ok = len(d) == len(main_) == f.drugs.iloc[-1] and (d.direction_group == "only consistent").sum() == len(only)
    rk = pd.read_csv(RK / "ranking_run_table.csv")
    mp = pd.read_csv(RUNS / "metapath_run_summary.csv").set_index("run")
    rank_ok = all(r.drugs_with_BB_routes - r.removed_liver_contraindicated == r.final_ranked_drugs
                  and r.drugs_with_BB_routes == mp.loc[r.run, "drugs_with_BB_routes"]
                  and r.final_drugs_treating_a_liver_disease == r.strict_liver_treating_drugs for r in rk.itertuples())
    for r in rk.itertuples():
        fin = pd.read_csv(RK / f"{r.run}_FINAL_ranked_drugs.csv")
        st = pd.read_csv(RK / f"{r.run}_strict_liver_treating_drugs.csv")
        rank_ok &= len(fin) == r.final_ranked_drugs and len(st) == r.strict_liver_treating_drugs
    before = (ROOT / "README.md").read_text(encoding="utf-8")
    subprocess.run([sys.executable, str(ROOT / "code/make_readme.py")], cwd=ROOT, check=True, capture_output=True)
    readme_ok = (ROOT / "README.md").read_text(encoding="utf-8") == before
    return {"funnel_steps_add_up": bool(steps_ok), "funnel_monotonic": bool(mono),
            "direction_groups_match_main_list": bool(dir_ok), "ranking_counts_match_files": bool(rank_ok),
            "readme_equals_fresh_build_from_files": readme_ok,
            "pass": bool(steps_ok and mono and dir_ok and rank_ok and readme_ok)}


def check5():
    files = [RUNS / f.format(r=r) for r in RUN_NAMES for f in ROUTE_FILES]
    files += [RK / f.format(r=r) for r in RUN_NAMES for f in RANK_FILES]
    files += [RK / "top25_comparison_FINAL.csv", RK / "top25_comparison_strict.csv"]
    before = {str(f.relative_to(ROOT)): table_hash(f) for f in files}
    backup = ROOT / "logs/before_rerun"
    if backup.exists():
        shutil.rmtree(backup)
    backup.mkdir(parents=True)
    for d in (RUNS, RK):
        shutil.move(str(d), str(backup / d.name))
        d.mkdir(parents=True)
    for sc in ["metapaths.py", "ranking.py"]:
        r = subprocess.run([sys.executable, str(ROOT / "code" / sc)], cwd=ROOT, capture_output=True, text=True)
        (LOGS / f"rerun_{sc.replace('.py', '')}.log").write_text(r.stdout + r.stderr, encoding="utf-8")
        if r.returncode:
            raise SystemExit(f"rerun of {sc} failed; previous outputs kept in {backup}")
    after = {k: table_hash(ROOT / k) for k in before}
    diff = [f for f in before if before[f] != after[f]]
    if not diff:
        shutil.rmtree(backup)
    return {"files_compared": len(before), "identical": not diff, "differing_files": diff,
            "note": "content compared after sorting rows; generation timestamps excluded", "pass": not diff}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res_path = OUT / "audit.json"
    res = json.loads(res_path.read_text(encoding="utf-8")) if res_path.exists() else {}
    if "--rerun" in sys.argv:
        res["5_rerun_identical"] = {"run": now(), **check5()}
    else:
        res.update({"1_excluded_names_in_code": {"run": now(), **check1()},
                    "2_graph_layers": {"run": now(), **check2()},
                    "3_gene_star_in_seed_file": {"run": now(), **check3()},
                    "4_numbers": {"run": now(), **check4()}})
    res_path.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
    lines = ["# Final audit", ""]
    for k, v in sorted(res.items()):
        lines.append(f"## {k}: {'PASS' if v['pass'] else 'FAIL'} ({v['run']})")
        if k.startswith("1"):
            lines.append(f"Searched {v['files_searched']} files for {len(v['patterns'])} patterns; hits: {len(v['hits'])}")
            lines += [f"- {h}" for h in v["hits"]]
        elif k.startswith("2"):
            lines.append(f"Database {v['database']}. TRACKS_* types: {v['tracks_relationship_types'] or 'none'}; "
                         f"expression-like labels: {v['expression_like_labels'] or 'none'}; "
                         f"Gene nodes with expression-like properties: {v['gene_nodes_with_expression_like_properties']}.")
            lines.append("Labels: " + ", ".join(f"{a} {b:,}" for a, b in v["labels"].items()))
            lines.append("Relationship types: " + ", ".join(f"{a} {b:,}" for a, b in v["relationship_types"].items()))
        elif k.startswith("3"):
            lines += [f"- {r}: {x['gene_star_in_routes']} Gene* in routes, outside the allowed set: {len(x['not_in_allowed_set'])}"
                      for r, x in v["runs"].items()]
        elif k.startswith("4"):
            lines += [f"- {a}: {b}" for a, b in v.items() if a not in ("run", "pass")]
        elif k.startswith("5"):
            lines.append(f"{v['files_compared']} files compared; identical: {v['identical']}; differing: {v['differing_files'] or 'none'}")
        lines.append("")
    (OUT / "audit_report.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
