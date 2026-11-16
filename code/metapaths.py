"""Metapaths B and BB for the only-consistent drugs, with the seed genes as Gene*.

Metapaths:
  B : Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
  BB: Drug -ACTS_ON-> Gene1 -INTERACTS_WITH- Gene2 -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
Rules:
  1. routes start only from the drug's direction-consistent targets (direction filter output)
  2. in-between genes have <= 500 INTERACTS_WITH partners (B: the target; BB: Gene1 and Gene2)
  3. pathway size 5-200 genes (IN_PATHWAY members)
  4. Gene1, Gene2 and Gene* are different genes, and the drug has no TREATS / OFF_LABEL_FOR link to
     NAFLD (EFO_0003095), NASH (EFO_1001249) or cirrhosis of liver (EFO_0001422)
BB is run only for drugs with kept B routes. Genes are matched by Ensembl ID (node_id).
Stage of a route = main_stage of its Gene* from the seed file (GWAS add-on genes: 0_genetic_risk).

Runs:
  RunA  all 280 seed genes (main result)
  RunB  228 histology genes (main_stage != 0_genetic_risk)
  RunC  genes with 3_fibrosis in all_stages
  RunD  280 seed genes + GWAS Catalog add-on genes (stage 0_genetic_risk, no direction)

Outputs: results/Metapath_runs/<run>_all_routes.csv (every B route + KEEP / drop reason), <run>.xlsx (kept B routes),
<run>_BB.csv and <run>_BB.xlsx (Summary + BB routes), metapath_run_summary.csv.

    python code/metapaths.py              # all runs
    python code/metapaths.py RunA RunB    # selected runs
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import NAFLD, ROOT, driver, now, read_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
GWAS = ROOT / "seed_genes/seed_genes_optional_GWAS_Catalog.csv"
DRUGS = ROOT / "results/Direction_consistent_drugs/only_consistent_drugs.csv"
LINKS = ROOT / "results/Direction_consistent_drugs/direction_consistent_links.csv"
OUT = ROOT / "results/Metapath_runs"
NAFLD_NAME = "non-alcoholic fatty liver disease (NAFLD)"
KNOWN_LIVER = ["EFO_0003095", "EFO_1001249", "EFO_0001422"]
HUB_MAX = 500
SIZE_MIN, SIZE_MAX = 5, 200

B_QUERY = """
UNWIND $pairs AS pr
MATCH (:Drug {node_id: pr.drug_id})-[:ACTS_ON]->(t:Gene {node_id: pr.gene_id})
MATCH (t)-[:INTERACTS_WITH]-(g:Gene) WHERE g.node_id IN $genes AND g <> t
WITH DISTINCT pr, t, g
MATCH (g)-[:IN_PATHWAY]->(p:Pathway)-[:IMPLICATED_IN]->(:Disease {node_id: $nafld})
RETURN pr.drug AS drug, pr.drug_id AS drug_id, t.node_id AS target_id, t.symbol AS target,
       COUNT { (t)-[:INTERACTS_WITH]-() } AS partner_n_partners,
       g.node_id AS gene_star_id, g.symbol AS gene_star, p.node_id AS pathway_id, p.name AS pathway,
       coalesce(p.source, 'Reactome') AS pathway_library, COUNT { (p)<-[:IN_PATHWAY]-(:Gene) } AS pathway_size
"""
PATHWAYS = """
MATCH (p:Pathway)-[:IMPLICATED_IN]->(:Disease {node_id: $nafld})
WITH p, COUNT { (p)<-[:IN_PATHWAY]-(:Gene) } AS size WHERE $smin <= size <= $smax
RETURN p.node_id AS id, p.name AS name, coalesce(p.source, 'Reactome') AS library, size
"""
BB_QUERY = """
UNWIND $pairs AS pr
MATCH (:Drug {node_id: pr.drug_id})-[:ACTS_ON]->(t:Gene {node_id: pr.gene_id})
WITH pr, t, COUNT { (t)-[:INTERACTS_WITH]-() } AS t_n WHERE t_n <= $hub
MATCH (t)-[:INTERACTS_WITH]-(h:Gene) WHERE h <> t
WITH DISTINCT pr, t, t_n, h
WITH pr, t, t_n, h, COUNT { (h)-[:INTERACTS_WITH]-() } AS h_n WHERE h_n <= $hub
MATCH (h)-[:INTERACTS_WITH]-(g:Gene) WHERE g.node_id IN $genes AND g <> t AND g <> h
WITH DISTINCT pr, t, t_n, h, h_n, g
MATCH (g)-[:IN_PATHWAY]->(p:Pathway) WHERE p.node_id IN $pathways
RETURN pr.drug AS drug, pr.drug_id AS drug_id, t.node_id AS gene1_id, t.symbol AS gene1, t_n AS gene1_n_partners,
       h.node_id AS gene2_id, h.symbol AS gene2, h_n AS gene2_n_partners, g.node_id AS gene_star_id, g.symbol AS gene_star,
       EXISTS { (t)-[:INTERACTS_WITH]-(g) } AS gene1_touches_gene_star, p.node_id AS pathway_id
"""


def gene_sets():
    s = pd.read_csv(SEEDS)
    s = s.rename(columns={"direction": "gene_star_direction"})[["gene", "ensembl_id", "main_stage", "all_stages", "gene_star_direction"]]
    g = pd.read_csv(GWAS)
    add = g[~g.ensembl_id.isin(s.ensembl_id)][["gene", "ensembl_id"]].assign(
        main_stage="0_genetic_risk", all_stages="0_genetic_risk", gene_star_direction="not_applicable")
    return {
        "RunA": ("RunA_all_seed_genes", s),
        "RunB": ("RunB_histology_genes", s[s.main_stage != "0_genetic_risk"]),
        # all_stages is joined with "; ", so strip each token before the comparison (an exact match on the
        # unstripped tokens only found genes whose FIRST stage is fibrosis: 122 instead of the 146 with fibrosis evidence)
        "RunC": ("RunC_fibrosis_genes", s[s.all_stages.str.split(";").map(lambda x: "3_fibrosis" in [t.strip() for t in x])]),
        "RunD": ("RunD_seed_plus_GWAS", pd.concat([s, add.assign(gene_star_source="GWAS Catalog add-on")], ignore_index=True)),
    }


def run_one(ses, name, genes, drugs, pairs, known):
    gidx = genes.set_index("ensembl_id")
    print(f"\n{name}: drugs {len(drugs)} | drug-target pairs {len(pairs)} | Gene* {len(genes)}", flush=True)
    rows = []
    for i in range(0, len(pairs), 50):
        rows += ses.run(B_QUERY, pairs=pairs[i:i + 50], genes=list(gidx.index), nafld=NAFLD).data()
    d = pd.DataFrame(rows, columns=["drug", "drug_id", "target_id", "target", "partner_n_partners", "gene_star_id",
                                    "gene_star", "pathway_id", "pathway", "pathway_library", "pathway_size"])
    d["run"] = name
    d["metapath"] = "B"
    d["direction_check"] = "direction-consistent"          # only consistent targets are used (rule 1)
    d["partner_gene"] = d.target
    d["gene_star_direction"] = d.gene_star_id.map(gidx.gene_star_direction)
    d["stage"] = d.gene_star_id.map(gidx.main_stage)
    d["route_id"] = ["|".join(map(str, t)) for t in zip(d.metapath, d.drug_id, d.target_id, d.gene_star_id, d.pathway_id)]
    why = []
    for r in d.itertuples():
        w = []
        if r.partner_n_partners > HUB_MAX:
            w.append(f"hub in-between gene (>{HUB_MAX} partners)")
        if not SIZE_MIN <= r.pathway_size <= SIZE_MAX:
            w.append(f"pathway size outside {SIZE_MIN}-{SIZE_MAX}")
        if known.get(r.drug_id, 0):
            w.append("drug already treats a liver disease (NAFLD/NASH/cirrhosis)")
        why.append("; ".join(w) or "KEEP")
    d["result"] = why
    d = d.drop_duplicates("route_id")
    d.sort_values(["result", "drug", "gene_star", "pathway"]).to_csv(OUT / f"{name}_all_routes.csv", index=False)
    kept = d[d.result == "KEEP"]
    cols = ["run", "route_id", "metapath", "drug", "drug_id", "target", "target_id", "direction_check", "partner_gene",
            "partner_n_partners", "gene_star", "gene_star_id", "gene_star_direction", "stage", "pathway", "pathway_id",
            "pathway_library", "pathway_size"]
    with pd.ExcelWriter(OUT / f"{name}.xlsx") as xw:
        (kept[cols] if len(kept) else pd.DataFrame([{"run": name, "route_id": "no kept routes"}], columns=cols)).sort_values(
            ["drug", "gene_star", "pathway"]).to_excel(xw, sheet_name="B", index=False)
    drop_reasons = d.result.str.split("; ").explode().value_counts().to_dict()
    res = {"run": name, "gene_star": len(genes), "drugs_in": len(drugs), "B_routes_found": len(d),
           "B_routes_kept": len(kept), "drugs_with_B_routes": int(kept.drug_id.nunique()), "B_drop_reasons": drop_reasons}
    print(f"  B: found {len(d)} | kept {len(kept)} | drugs with kept B routes {kept.drug_id.nunique()} | {drop_reasons}", flush=True)
    if kept.empty:
        res.update(BB_routes=0, drugs_with_BB_routes=0)
        return res

    b_drugs = set(kept.drug_id)
    b_routes = kept.groupby("drug_id").size()
    bpairs = [p for p in pairs if p["drug_id"] in b_drugs]
    pw = pd.DataFrame(ses.run(PATHWAYS, nafld=NAFLD, smin=SIZE_MIN, smax=SIZE_MAX).data()).set_index("id")
    rows = []
    for i in range(0, len(bpairs), 10):
        rows += ses.run(BB_QUERY, pairs=bpairs[i:i + 10], genes=list(gidx.index), hub=HUB_MAX, pathways=list(pw.index)).data()
    r = pd.DataFrame(rows).drop_duplicates(["drug_id", "gene1_id", "gene2_id", "gene_star_id", "pathway_id"]) if rows else pd.DataFrame()
    res.update(BB_routes=len(r), drugs_with_BB_routes=int(r.drug_id.nunique()) if len(r) else 0, BB_pathways_eligible=len(pw))
    print(f"  BB: routes {len(r)} | drugs with BB routes {res['drugs_with_BB_routes']} of {len(b_drugs)}", flush=True)
    if r.empty:
        return res
    r["pathway"] = r.pathway_id.map(pw.name)
    r["pathway_library"] = r.pathway_id.map(pw.library)
    r["pathway_size"] = r.pathway_id.map(pw["size"])
    r["gene_star_direction"] = r.gene_star_id.map(gidx.gene_star_direction)
    r["stage"] = r.gene_star_id.map(gidx.main_stage)
    r["disease"] = NAFLD_NAME
    r["also_in_B"] = [f"B|{a}|{b}|{c}|{p}" in set(kept.route_id) and b == g1
                      for a, b, c, p, g1 in zip(r.drug_id, r.gene1_id, r.gene_star_id, r.pathway_id, r.gene1_id)]
    r = r.sort_values(["drug", "gene1", "gene2", "gene_star", "pathway"])
    r.to_csv(OUT / f"{name}_BB.csv", index=False)
    s = r.groupby(["drug_id", "drug"]).size().rename("BB routes").reset_index()
    s["B routes"] = s.drug_id.map(b_routes).fillna(0).astype(int)
    no = sorted(b_drugs - set(r.drug_id))
    if no:
        nm = kept.drop_duplicates("drug_id").set_index("drug_id").drug
        s = pd.concat([s, pd.DataFrame({"drug_id": no, "drug": [nm[x] for x in no], "BB routes": 0,
                                        "B routes": [int(b_routes[x]) for x in no]})])
    s = s.sort_values(["BB routes", "drug"], ascending=[False, True]).reset_index(drop=True)
    s.insert(0, "No.", s.index + 1)
    s.insert(1, "Drug -> Disease", s.drug + " -> NAFLD")
    s["Disease"] = NAFLD_NAME
    s = s.rename(columns={"drug": "Drug", "drug_id": "Drug_ID"})
    simple = r[["drug", "gene1", "gene2", "gene_star", "stage", "pathway", "disease"]]
    simple.columns = ["Drug", "Gene 1", "Gene 2", "Gene*", "Gene* stage", "Pathway", "Disease"]
    with pd.ExcelWriter(OUT / f"{name}_BB.xlsx") as xw:
        s.to_excel(xw, sheet_name="Summary", index=False)
        step = 1_000_000
        for k, i in enumerate(range(0, max(len(simple), 1), step), 1):
            simple.iloc[i:i + step].to_excel(xw, sheet_name="BB" if k == 1 else f"BB ({k})", index=False)
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    drugs = pd.read_csv(DRUGS)
    links = pd.read_csv(LINKS)
    links = links[links.Drug_ID.isin(drugs.Drug_ID)]
    pairs = (links[["Drug", "Drug_ID", "Gene_ID"]].drop_duplicates()
             .rename(columns={"Drug": "drug", "Drug_ID": "drug_id", "Gene_ID": "gene_id"}).to_dict("records"))
    assert links.direction_check.eq("direction-consistent").all()
    sets = gene_sets()
    todo = sys.argv[1:] or list(sets)
    drv = driver()
    results = []
    with read_session(drv) as ses:
        if not ses.run("MATCH ()-[r:IMPLICATED_IN]->(:Disease {node_id: $d}) RETURN count(r) > 0", d=NAFLD).single()[0]:
            raise SystemExit("no Pathway -IMPLICATED_IN-> NAFLD links (step 0.7)")
        known = {r["id"]: r["n"] for r in ses.run(
            """UNWIND $ids AS id MATCH (d:Drug {node_id: id})
               OPTIONAL MATCH (d)-[:TREATS|OFF_LABEL_FOR]->(x:Disease) WHERE x.node_id IN $liver
               RETURN id, count(x) AS n""", ids=list(drugs.Drug_ID), liver=KNOWN_LIVER)}
        for key in todo:
            name, genes = sets[key]
            results.append({"generated": now(), **run_one(ses, name, genes, drugs, pairs, known)})
    drv.close()
    summ = OUT / "metapath_run_summary.csv"
    old = pd.read_csv(summ) if summ.exists() else pd.DataFrame()
    new = pd.DataFrame(results)
    if len(old):
        old = old[~old.run.isin(new.run)]
    pd.concat([old, new], ignore_index=True).sort_values("run").to_csv(summ, index=False)
    print("\n" + pd.read_csv(summ).drop(columns=["generated"]).to_string(index=False))


if __name__ == "__main__":
    main()
