"""Download liver-disease GWAS hits from the GWAS Catalog REST API v2.

Keeps genome-wide significant associations (p <= 5e-8) and their mapped genes.
Output columns: efo, trait, group, gene, p, pubmed, snp

Usage:
    python scripts/fetch_gwas_catalog.py source/gwas_catalog_liver_hits.csv
"""
import json, re, sys, urllib.request

import pandas as pd

API = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations?efo_id={}&size=200&page={}"
TRAITS = {
    "MONDO_0004790": ("fatty liver disease", "core"),
    "MONDO_0013209": ("MASLD", "core"),
    "EFO_0008421": ("NAFLD severity", "core"),
    "MONDO_0007027": ("MASH", "core"),
    "MONDO_0005155": ("cirrhosis of liver", "core"),
    "EFO_0010576": ("liver fibrosis measurement", "core"),
    "EFO_0010821": ("liver fat measurement", "extended"),
}

rows = []
for efo, (name, grp) in TRAITS.items():
    page = 0
    while True:
        with urllib.request.urlopen(API.format(efo, page), timeout=120) as r:
            d = json.load(r)
        for a in d.get("_embedded", {}).get("associations", []):
            for g in a.get("mapped_genes") or []:
                for gg in re.split(r"[,;]\s*", g):
                    if gg and gg != "-":
                        rows.append(dict(efo=efo, trait=name, group=grp, gene=gg.strip(),
                                         p=a.get("p_value"), pubmed=a.get("pubmed_id"),
                                         snp=";".join(x["rs_id"] for x in a.get("snp_allele") or [])))
        page += 1
        if page >= d["page"]["totalPages"]:
            break
    print(efo, name, d["page"]["totalElements"], "associations")

out = pd.DataFrame(rows)
out = out[out.p <= 5e-8]
out.to_csv(sys.argv[1], index=False)
print(len(out), "significant gene-trait rows,", out.gene.nunique(), "genes")
