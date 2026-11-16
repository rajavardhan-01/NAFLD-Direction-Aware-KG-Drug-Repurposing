# Seed genes (Gene*)

The seed genes are the disease-side anchor of every metapath route (Gene\*). All come from public, citable sources.

## Main file: `seed_genes_main_Ryaboshapkina2017.csv` (280 genes)

Source: Ryaboshapkina M, Hammar M. *Human hepatic gene expression signature of non-alcoholic fatty liver disease
progression, a meta-analysis.* Sci Rep 2017;7:12361. doi:10.1038/s41598-017-10930-w (open access, CC BY 4.0),
Supplementary Table S2.

- 218 genes whose liver expression tracks NAFLD histology, each replicated with the same direction in at least three
  independent public studies.
- 62 genes with a reported genetic association with NAFLD, compiled by the same authors (52 of them have no histology
  signal and are labelled `0_genetic_risk`).

| Column | Meaning |
|---|---|
| `gene` | current HGNC symbol (6 symbols from 2016 updated, e.g. PAPD7 -> TENT4A) |
| `ensembl_id` | Ensembl gene ID (= OptimusKG gene ID) |
| `main_stage` | earliest disease stage the gene is linked to |
| `all_stages` | every stage the gene is linked to |
| `evidence` | histology traits the paper reported for the gene |
| `direction` | `increasing` / `decreasing` as disease gets worse; `not_applicable` for genetic-risk genes |
| `n_gene_partners`, `n_pathways`, `n_drugs_targeting` | connections in OptimusKG |
| `usable_for_metapath` | at least one gene partner and one pathway in OptimusKG |
| `also_in_GWAS_Catalog` | gene also among current GWAS Catalog liver hits |

### Disease stages

| Stage | Trait in the paper | Genes (main stage) | Increasing / decreasing |
|---|---|---|---|
| `0_genetic_risk` | genetic association | 52 | n/a |
| `1_steatosis` | degree of steatosis | 13 | 10 / 3 |
| `2_NASH` | NASH vs simple steatosis, inflammation, NAS | 83 | 56 / 27 |
| `3_fibrosis` | fibrosis progression | 122 | 93 / 29 |
| `4_HCC` | progression to HCC | 10 | 1 / 9 |

## Optional add-on: `seed_genes_optional_GWAS_Catalog.csv` (310 genes)

Genome-wide significant (p <= 5e-8) genes for fatty liver disease, MASLD, MASH, NAFLD severity, cirrhosis and liver
fibrosis from the GWAS Catalog (downloaded 2026-11-11); pseudogenes, antisense and lncRNA names removed. Adds newer
risk genes (e.g. PNPLA3, TM6SF2, HSD17B13, MTARC1, GCKR). They have no stage or direction and are used as
`0_genetic_risk` in RunD only.

`seed_genes_main_long_format.csv` has one row per gene x trait.

## Rebuild

```
python scripts/fetch_gwas_catalog.py source/gwas_catalog_liver_hits_2026-11-11.csv
python scripts/build_seeds.py source/Ryaboshapkina2017_SciRep_SupplementaryTables.xls source/gwas_catalog_liver_hits_2026-11-11.csv <optimuskg-parquet-folder> .
```

- `source/Ryaboshapkina2017_SciRep_SupplementaryTables.xls`: original supplement (CC BY 4.0)
- `source/gwas_catalog_liver_hits_2026-11-11.csv`: raw GWAS Catalog download
