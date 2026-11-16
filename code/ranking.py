"""Ranking: the final ranked drug list for each run, plus a strict liver-treating view, stage tables and
the top-25 comparison across runs.

Input: each run's BB routes (results/Metapath_runs/<run>_BB.csv/.xlsx). Output folder: results/Ranking/.

FINAL RESULT (per run)
  <run>_FINAL_ranked_drugs.csv/.xlsx    every drug with >= 1 BB route to NAFLD, minus drugs CONTRAINDICATED_IN any
                                        liver disease, ranked by BB routes (ties share the lowest rank)
  <run>_FINAL_top_drugs_stages.xlsx     ranks 1-25 (ties kept): BB routes per Gene* stage 0_genetic_risk ... 4_HCC
  top25_comparison_FINAL.csv/.md        top 25 of the final lists across runs

STRICT VIEW (sensitivity check; keeps only drugs that already treat a liver disease)
  <run>_strict_liver_treating_drugs.csv/.xlsx   final-list drugs with >= 1 TREATS link to a liver disease
                                                (Task 4 liver set); off-label use alone does not count
  <run>_strict_top_drugs_stages.xlsx, top25_comparison_strict.csv/.md

Also: <run>_direct_diseases.xlsx (every TREATS / OFF_LABEL_FOR / CONTRAINDICATED_IN disease of the BB drugs) and
ranking_run_table.csv. Read-only Neo4j.

    python code/ranking.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, driver, now, read_session  # noqa: E402

RUNS_DIR = ROOT / "results/Metapath_runs"
OUT = ROOT / "results/Ranking"
LIVER_SET = ROOT / "results/Task4/liver_disease_set.csv"
RUNS = {"RunA": "RunA_all_seed_genes", "RunB": "RunB_histology_genes", "RunC": "RunC_fibrosis_genes",
        "RunD": "RunD_seed_plus_GWAS"}
STAGES = ["0_genetic_risk", "1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]
TOP_RANK = 25
LINK = {"TREATS": "treats", "OFF_LABEL_FOR": "off-label for", "CONTRAINDICATED_IN": "contraindicated in"}


def join(s):
    return "; ".join(sorted(set(s.dropna())))


def rank(df):
    df = df.sort_values(["BB routes", "Drug"], ascending=[False, True]).copy()
    df.insert(0, "Rank", df["BB routes"].rank(method="min", ascending=False).astype(int))
    return df.rename(columns={"BB routes": "BB routes to NAFLD"})


def rank_run(ses, run, liver_ids):
    summ = pd.read_excel(RUNS_DIR / f"{run}_BB.xlsx", sheet_name="Summary")
    routes = pd.read_csv(RUNS_DIR / f"{run}_BB.csv")
    bb = summ[summ["BB routes"] > 0][["Drug_ID", "Drug", "BB routes", "B routes"]]
    links = pd.DataFrame(ses.run("""UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:TREATS|OFF_LABEL_FOR|CONTRAINDICATED_IN]->(x:Disease)
        RETURN id AS Drug_ID, type(r) AS rel, x.node_id AS disease_id, coalesce(x.name, x.node_id) AS Disease""",
                                 ids=list(bb.Drug_ID)).data(), columns=["Drug_ID", "rel", "disease_id", "Disease"])
    links = links.drop_duplicates(["Drug_ID", "rel", "disease_id"]).merge(bb[["Drug_ID", "Drug"]], on="Drug_ID")
    links["Link"] = links.rel.map(LINK)
    links["is_liver"] = links.disease_id.isin(liver_ids)
    lv = links[links.is_liver]

    S = bb.rename(columns={"BB routes": "Indirect routes to NAFLD"}).copy()
    for rel, col in [("TREATS", "Treats"), ("OFF_LABEL_FOR", "Off-label for"), ("CONTRAINDICATED_IN", "Contraindicated in")]:
        sub = links[links.rel == rel].groupby("Drug_ID").Disease
        S[f"{col.split()[0]} (n)"] = S.Drug_ID.map(sub.nunique()).fillna(0).astype(int)
        S[col] = S.Drug_ID.map(sub.agg(join)).fillna("")
    S["Liver-related direct links"] = S.Drug_ID.map(lv.assign(t=lv.Link + " " + lv.Disease).groupby("Drug_ID").t.agg(join)).fillna("")
    S = S.sort_values(["Indirect routes to NAFLD", "Drug"], ascending=[False, True]).reset_index(drop=True)
    S.insert(0, "No.", S.index + 1)
    with pd.ExcelWriter(OUT / f"{run}_direct_diseases.xlsx") as xw:
        S[["No.", "Drug", "Drug_ID", "Indirect routes to NAFLD", "Treats (n)", "Off-label (n)", "Contraindicated (n)",
           "Treats", "Off-label for", "Contraindicated in", "Liver-related direct links", "B routes"]].to_excel(
            xw, sheet_name="Summary", index=False)
        links.sort_values(["Drug", "Link", "Disease"])[["Drug", "Link", "Disease", "disease_id"]].to_excel(
            xw, sheet_name="Direct links", index=False)
        lv.sort_values(["Drug", "Link", "Disease"])[["Drug", "Link", "Disease", "disease_id"]].to_excel(
            xw, sheet_name="Liver links", index=False)

    treats_liver = lv[lv.rel == "TREATS"].groupby("Drug_ID").Disease.agg(join)
    ci_liver = lv[lv.rel == "CONTRAINDICATED_IN"].groupby("Drug_ID").Disease.agg(join)
    removed = bb[bb.Drug_ID.isin(ci_liver.index)].copy()
    removed["Liver diseases it treats"] = removed.Drug_ID.map(treats_liver)
    removed["Liver contraindications"] = removed.Drug_ID.map(ci_liver)
    removed = removed.sort_values("BB routes", ascending=False).rename(columns={"BB routes": "BB routes to NAFLD"})

    # FINAL list: all BB drugs minus liver-contraindicated drugs
    F = bb[~bb.Drug_ID.isin(ci_liver.index)].copy()
    F["Treats a liver disease"] = F.Drug_ID.isin(treats_liver.index)
    F["Liver diseases it treats"] = F.Drug_ID.map(treats_liver)
    F["Current indications (TREATS)"] = F.Drug_ID.map(links[links.rel == "TREATS"].groupby("Drug_ID").Disease.agg(join))
    F = rank(F)[["Rank", "Drug", "Drug_ID", "BB routes to NAFLD", "B routes", "Treats a liver disease",
                 "Liver diseases it treats", "Current indications (TREATS)"]]
    F.to_csv(OUT / f"{run}_FINAL_ranked_drugs.csv", index=False)
    with pd.ExcelWriter(OUT / f"{run}_FINAL_ranked_drugs.xlsx") as xw:
        F.to_excel(xw, sheet_name="Final ranked drugs", index=False)
        removed.to_excel(xw, sheet_name="Removed (liver CI)", index=False)
    stage_tables(F, routes, OUT / f"{run}_FINAL_top_drugs_stages.xlsx")

    # STRICT view: final-list drugs that already treat a liver disease, re-ranked
    T = rank(bb[bb.Drug_ID.isin(treats_liver.index) & ~bb.Drug_ID.isin(ci_liver.index)].assign(
        **{"Liver diseases it treats": lambda d: d.Drug_ID.map(treats_liver)}))
    T = T[["Rank", "Drug", "Drug_ID", "BB routes to NAFLD", "B routes", "Liver diseases it treats"]]
    T.to_csv(OUT / f"{run}_strict_liver_treating_drugs.csv", index=False)
    with pd.ExcelWriter(OUT / f"{run}_strict_liver_treating_drugs.xlsx") as xw:
        T.to_excel(xw, sheet_name="Strict liver-treating", index=False)
        removed[removed["Liver diseases it treats"].notna()].to_excel(xw, sheet_name="Removed (liver CI)", index=False)
    if len(T):
        stage_tables(T, routes, OUT / f"{run}_strict_top_drugs_stages.xlsx")

    info = {"run": run, "drugs_with_B_routes": len(summ), "drugs_with_BB_routes": len(bb),
            "removed_liver_contraindicated": len(removed), "final_ranked_drugs": len(F),
            "final_rank_le_25": int((F.Rank <= TOP_RANK).sum()),
            "final_drugs_treating_a_liver_disease": int(F["Treats a liver disease"].sum()),
            "strict_liver_treating_drugs": len(T),
            "strict_treat_liver_but_removed_ci": int(removed["Liver diseases it treats"].notna().sum())}
    (OUT / f"{run}_ranking.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    return info, F, T


def stage_tables(final, routes, path):
    top = final[final.Rank <= TOP_RANK]
    r = routes[routes.drug_id.isin(top.Drug_ID)].copy()
    st = r.groupby(["drug_id", "stage"]).size().unstack(fill_value=0).reindex(columns=STAGES, fill_value=0)
    D = top[["Rank", "Drug", "Drug_ID", "BB routes to NAFLD"]].rename(columns={"BB routes to NAFLD": "BB routes"}).set_index("Drug_ID")
    D["Gene* (n)"] = r.groupby("drug_id").gene_star_id.nunique()
    for s_ in STAGES:
        D[f"{s_} routes"] = st[s_]
    D["Main stage (most routes)"] = st.idxmax(axis=1)
    D["% routes in main stage"] = (100 * st.max(axis=1) / st.sum(axis=1)).round(1)
    D["All stages reached"] = st.apply(lambda row: ", ".join(c for c in STAGES if row[c] > 0), axis=1)
    D["Liver disease it treats"] = top.set_index("Drug_ID")["Liver diseases it treats"]
    D = D.reset_index()[["Rank", "Drug", "Drug_ID", "BB routes", "Gene* (n)"] + [f"{s_} routes" for s_ in STAGES]
                        + ["Main stage (most routes)", "% routes in main stage", "All stages reached", "Liver disease it treats"]]
    gs = r.groupby(["drug_id", "drug", "gene_star", "stage", "gene_star_direction"]).size().rename("BB routes").reset_index()
    gs = gs.rename(columns={"drug": "Drug", "gene_star": "Gene*", "stage": "Gene* stage", "gene_star_direction": "Direction"})
    gs["_rank"] = gs.drug_id.map(top.set_index("Drug_ID").Rank)
    gs = gs.sort_values(["_rank", "Drug", "BB routes"], ascending=[True, True, False]).drop(columns=["_rank", "drug_id"])
    tot = pd.DataFrame({"Stage": STAGES, "Drugs reaching it": [int((st[s_] > 0).sum()) for s_ in STAGES],
                        "Routes": [int(st[s_].sum()) for s_ in STAGES],
                        "Main stage for (n drugs)": [int((D["Main stage (most routes)"] == s_).sum()) for s_ in STAGES]})
    with pd.ExcelWriter(path) as xw:
        D.to_excel(xw, sheet_name="Drug stages", index=False)
        gs.to_excel(xw, sheet_name="Drug - Gene star detail", index=False)
        tot.to_excel(xw, sheet_name="Stage totals", index=False)


def compare(lists, stem, title):
    top = {k: f[f.Rank <= TOP_RANK].set_index("Drug_ID").Rank for k, f in lists.items()}
    liver_tx = {d: t for f in lists.values() for d, t in zip(f.Drug_ID, f["Liver diseases it treats"])}
    names = {d: n for f in lists.values() for d, n in zip(f.Drug_ID, f.Drug)}
    t = pd.DataFrame(top).reindex(columns=list(lists))
    t["in_all_runs"] = t[list(top)].notna().all(axis=1)
    t["in_n_runs"] = t[list(top)].notna().sum(axis=1)
    t["liver_diseases_it_treats"] = t.index.map(liver_tx)
    t = t.sort_values(["in_n_runs"] + list(top), ascending=[False] + [True] * len(top))
    t.index.name = "Drug_ID"
    t.insert(0, "Drug", t.index.map(names))
    t.to_csv(OUT / f"{stem}.csv")
    sets = {k: set(v.index) for k, v in top.items()}
    keys = list(sets)
    fmt = lambda x: "" if pd.isna(x) else str(int(x))  # noqa: E731
    over = [("All runs", len(set.intersection(*sets.values())))] + [
        (f"{a} & {b}", len(sets[a] & sets[b])) for i, a in enumerate(keys) for b in keys[i + 1:]]
    md = f"""# {title}

Top 25 = rank <= 25 by BB routes to NAFLD (ties included): {", ".join(f"{k} {len(v)}" for k, v in sets.items())} drugs.
Generated {now()}.

| Overlap | Drugs |
|---|---|
""" + "\n".join(f"| {a} | {n} |" for a, n in over) + f"""

| Drug | Drug_ID | {" | ".join("Rank " + k for k in keys)} | In all runs | Liver diseases it treats |
|---|---|{"---|" * len(keys)}---|---|
""" + "\n".join(f"| {r.Drug} | {d} | {' | '.join(fmt(r[k]) for k in keys)} | {'yes' if r.in_all_runs else ''} | "
                f"{'' if pd.isna(r.liver_diseases_it_treats) else r.liver_diseases_it_treats} |" for d, r in t.iterrows()) + "\n"
    (OUT / f"{stem}.md").write_text(md, encoding="utf-8")
    return md


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    liver = pd.read_csv(LIVER_SET)
    liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
    drv = driver()
    infos, finals, stricts = [], {}, {}
    with read_session(drv) as ses:
        for key, run in RUNS.items():
            info, F, T = rank_run(ses, run, liver_ids)
            infos.append(info)
            finals[key], stricts[key] = F, T
            print(json.dumps(info))
    drv.close()
    tab = pd.DataFrame(infos)
    tab.to_csv(OUT / "ranking_run_table.csv", index=False)
    md = compare(finals, "top25_comparison_FINAL",
                 "FINAL result: top 25 ranked drugs across runs (all drugs with BB routes; liver-contraindicated drugs removed)")
    compare(stricts, "top25_comparison_strict",
            "Strict view: top 25 drugs that already treat a liver disease, across runs")
    print(tab.to_string(index=False))
    print(md)


if __name__ == "__main__":
    main()
