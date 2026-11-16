"""Build public seed genes from Ryaboshapkina & Hammar 2017 (Sci Rep, Table S2)
plus current GWAS Catalog liver hits, and check connectivity in OptimusKG."""
import re, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

XLS, GWAS, OKG, OUTDIR = sys.argv[1:5]

d = pd.read_excel(XLS, sheet_name="Metaanalysis_NAFLD", header=2)
# 2016 HGNC symbols renamed since; keep the paper's name in a separate column
RENAMED = {"FAM129B": "NIBAN2", "FAM213B": "PRXL2B", "GUCY1A3": "GUCY1A1",
           "H2AFY": "MACROH2A1", "HIST1H2AC": "H2AC6", "PAPD7": "TENT4A"}
d["symbol_in_paper"] = d.GeneSymbol
d["GeneSymbol"] = d.GeneSymbol.replace(RENAMED)

# disease stages in progression order, with the three studies the paper used for each
STAGES = [
    ("1_steatosis", "highConf_degree_of_SS",
     ["GSE48452_degree_of_steatosis", "GSE61260_NAFLD_degree_of_steatosis", "GSE89632_degree_of_steatosis"]),
    ("2_NASH", "highConf_NASH_vs_SS",
     ["GSE33814_nash_vs_simplesteatosis", "GSE48452_nash_vs_simplesteatosis", "GSE89632_nash_vs_simplesteatosis"]),
    ("2_NASH", "highConf_inflammation",
     ["GSE48452_inflammation", "GSE61260_NAFLD_inflammation", "GSE89632_inflammation"]),
    ("2_NASH", "highConf_NAS",
     ["GSE48452_nas", "GSE61260_NAFLD_nas", "GSE89632_nas"]),
    ("3_fibrosis", "highConf_fibrosis_in_NAFLD",
     ["GSE48452_fibrosis_progression", "GSE49541_fibrosis_progression", "GSE89632_fibrosis_progression"]),
    ("4_HCC", "highConf_progression_to_HCC",
     ["GSE14323_progression_to_HCC", "GSE54238_progression_to_HCC", "GSE6764_progression_to_HCC"]),
]

def sign(row, cols):
    v = pd.to_numeric(row[cols], errors="coerce").dropna()
    s = v.sum()
    return "increasing" if s > 0 else "decreasing" if s < 0 else "inconsistent"

rows = []
for _, r in d.iterrows():
    genetic = str(r.GeneCategory).startswith("genetic")
    hits = [(stg, flag, sign(r, cols)) for stg, flag, cols in STAGES if r.get(flag) == 1]
    if hits:
        for stg, flag, dirn in hits:
            rows.append(dict(gene=r.GeneSymbol, entrez=r.entrezids, stage=stg,
                             evidence=flag.replace("highConf_", ""), direction=dirn,
                             source="Ryaboshapkina2017_TableS2_histology"))
    elif genetic:
        rows.append(dict(gene=r.GeneSymbol, entrez=r.entrezids, stage="0_genetic_risk",
                         evidence="genetic association with NAFLD", direction="not_applicable",
                         source="Ryaboshapkina2017_TableS2_genetic"))
long = pd.DataFrame(rows)

# one row per gene: earliest stage it appears at = its "stage"; keep all stages listed
def collapse(g):
    g = g.sort_values("stage")
    first = g.iloc[0]
    dirs = set(g.direction) - {"not_applicable"}
    return pd.Series(dict(entrez=first.entrez, main_stage=first.stage,
                          all_stages="; ".join(dict.fromkeys(g.stage)),
                          evidence="; ".join(dict.fromkeys(g.evidence)),
                          direction=(dirs.pop() if len(dirs) == 1 else
                                     "mixed" if dirs else "not_applicable"),
                          source=first.source))
seeds = long.groupby("gene").apply(collapse, include_groups=False).reset_index()

# GWAS Catalog (current) genes as optional add-on
junk = re.compile(r"(^RNU|^RNA5S|^RN7S|P\d+$|-AS\d*$|-DT$|^LINC|^MIR|^SNOR|^Y_RNA|^LOC|-IT\d$|^ENSG)")
gw = pd.read_csv(GWAS)
gw = gw[(gw.group == "core") & ~gw.gene.str.contains(junk)]
gwg = gw.groupby("gene").agg(gwas_traits=("trait", lambda s: "; ".join(sorted(set(s)))),
                             gwas_best_p=("p", "min"),
                             gwas_pubmed=("pubmed", lambda s: "; ".join(sorted(set(map(str, s)))))).reset_index()

# OptimusKG connectivity
g = pd.read_parquet(f"{OKG}/gene.parquet", columns=["id", "properties"])
g["sym"] = g.properties.map(lambda p: p.get("symbol") or p.get("name") or "")
sym2id = dict(zip(g.sym, g.id))
gg = pq.read_table(f"{OKG}/gene_gene.parquet", columns=["from", "to"]).to_pandas()
deg = pd.concat([gg["from"], gg["to"]]).value_counts()
pg = pq.read_table(f"{OKG}/pathway_gene.parquet", columns=["from", "to"]).to_pandas()
pcol = "to" if pg["to"].astype(str).str.startswith("ENSG").mean() > 0.5 else "from"
pdeg = pg[pcol].value_counts()
dg = pq.read_table(f"{OKG}/drug_gene.parquet", columns=["from", "to"]).to_pandas()
dcol = "to" if dg["to"].astype(str).str.startswith("ENSG").mean() > 0.5 else "from"
ddeg = dg[dcol].value_counts()

def okg(df):
    df["ensembl_id"] = df.gene.map(sym2id)
    df["in_OptimusKG"] = df.ensembl_id.notna()
    df["n_gene_partners"] = df.ensembl_id.map(deg).fillna(0).astype(int)
    df["n_pathways"] = df.ensembl_id.map(pdeg).fillna(0).astype(int)
    df["n_drugs_targeting"] = df.ensembl_id.map(ddeg).fillna(0).astype(int)
    df["usable_for_metapath"] = (df.n_gene_partners > 0) & (df.n_pathways > 0)
    return df

seeds = okg(seeds)
seeds["also_in_GWAS_Catalog"] = seeds.gene.isin(gwg.gene)
gwg = okg(gwg)
gwg["also_in_Ryaboshapkina2017"] = gwg.gene.isin(seeds.gene)

seeds.sort_values(["main_stage", "gene"]).to_csv(f"{OUTDIR}/seed_genes_main_Ryaboshapkina2017.csv", index=False)
long.to_csv(f"{OUTDIR}/seed_genes_main_long_format.csv", index=False)
gwg.sort_values("gwas_best_p").to_csv(f"{OUTDIR}/seed_genes_optional_GWAS_Catalog.csv", index=False)

print("pathway gene column:", pcol, " drug gene column:", dcol)
print("MAIN seeds:", len(seeds))
print(seeds.groupby("main_stage").agg(n=("gene", "size"), in_okg=("in_OptimusKG", "sum"),
                                      usable=("usable_for_metapath", "sum")).to_string())
print(pd.crosstab(seeds.main_stage, seeds.direction).to_string())
print("missing from OptimusKG:", sorted(seeds.loc[~seeds.in_OptimusKG, "gene"]))
print("histology genes total:", (seeds.source.str.contains("histology")).sum(),
      " usable:", seeds.loc[seeds.source.str.contains("histology"), "usable_for_metapath"].sum())
print("genetic:", seeds.loc[seeds.main_stage == "0_genetic_risk", "gene"].tolist())
print("\nGWAS optional:", len(gwg), " in OKG:", gwg.in_OptimusKG.sum(), " usable:", gwg.usable_for_metapath.sum(),
      " already in main:", gwg.also_in_Ryaboshapkina2017.sum())
print("main fibrosis-stage genes:", sorted(seeds.loc[seeds.all_stages.str.contains("3_fibrosis"), "gene"])[:200])
