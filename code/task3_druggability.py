"""Task 3: druggability of the seed genes, overall and by stage.

Background: all genes in this graph with >= 1 INTERACTS_WITH partner (OptimusKG gene-gene network).
Seed genes without a partner are outside the background and are left out of the tests (listed).

Druggable (main definition) = pharmacological target: >= 1 ACTS_ON link that is
  an OptimusKG action other than ENZYME / TRANSPORTER / CARRIER, or a curated link added in step 0.4,
  or an OptimusKG link carrying a GtoPdb/DGIdb/ChEMBL direction from step 0.5.
Secondary definition: any drug (any ACTS_ON link).

Groups: all seed genes; histology genes (main_stage != 0_genetic_risk); each main_stage;
each histology stage x direction (increasing / decreasing).
Statistics: Fisher's exact test vs the rest of the background (odds ratio, 95% CI), BH-FDR over
groups; pairwise Fisher tests between the four histology stages, BH-FDR; 1,000 degree-matched permutations
(INTERACTS_WITH degree: degree-0 bin is empty here, then deciles), seed 42.

Outputs: results/Task3/.

    python code/task3_druggability.py
"""
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import false_discovery_control, fisher_exact  # noqa: E402
from scipy.stats.contingency import odds_ratio  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import ROOT, driver, now, read_session  # noqa: E402

SEEDS = ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv"
EDGES = ROOT / "results/Task2/drug_gene_edges.csv"
OUT = ROOT / "results/Task3"
N_PERM, SEED, N_DEGREE_BINS = 1000, 42, 10
PK_ACTIONS = ["ENZYME", "TRANSPORTER", "CARRIER"]
ADDED_LAYER = "added curated drug-target links"
STAGES = ["0_genetic_risk", "1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]
DEFS = {"pharm": ("pharm_connected", "pharmacological target (main)", "_pharmacological"),
        "any": ("drug_connected", "any drug link", "_any")}

QUERY = """
MATCH (g:Gene) WHERE (g)-[:INTERACTS_WITH]-(:Gene)
CALL (g) { OPTIONAL MATCH (g)-[:INTERACTS_WITH]-(h:Gene) WHERE h <> g RETURN count(DISTINCT h) AS ppi_degree }
CALL (g) { OPTIONAL MATCH (d:Drug)-[r:ACTS_ON]->(g)
           RETURN count(DISTINCT d) AS n_drugs,
                  count(DISTINCT CASE WHEN r.layer = $layer OR r.external_action_direction IS NOT NULL
                                           OR (r.action IS NOT NULL AND NOT r.action IN $pk) THEN d END) AS n_pharm_drugs }
RETURN g.ensembl_id AS ensembl_id, ppi_degree, n_drugs, n_pharm_drugs
"""


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = pd.read_csv(SEEDS)
    drv = driver()
    with read_session(drv) as s:
        bg = pd.DataFrame(s.run(QUERY, layer=ADDED_LAYER, pk=PK_ACTIONS).data())
    drv.close()
    bg = bg[bg.ppi_degree > 0].drop_duplicates("ensembl_id").reset_index(drop=True)
    bg["drug_connected"] = bg.n_drugs > 0
    bg["pharm_connected"] = bg.n_pharm_drugs > 0
    bg["degree_bin"] = pd.qcut(bg.ppi_degree.rank(method="first"), N_DEGREE_BINS, labels=False) + 1
    bg.to_csv(OUT / "background_genes_druggability.csv", index=False)

    outside = seeds[~seeds.ensembl_id.isin(bg.ensembl_id)]
    outside[["gene", "ensembl_id", "main_stage", "direction"]].to_csv(OUT / "seed_genes_outside_background.csv", index=False)
    s_in = seeds[seeds.ensembl_id.isin(bg.ensembl_id)]

    # consistency with Task 2
    edges = pd.read_csv(EDGES, low_memory=False)
    edges["pharmacological"] = ((edges.Link_origin != "OptimusKG") | edges.external_action_direction.notna()
                                | (~edges.Action_raw.isin(PK_ACTIONS) & (edges.Action_raw != "unknown")))
    bgi = bg.set_index("ensembl_id")
    in_seed = bg.ensembl_id.isin(s_in.ensembl_id)
    assert set(bg.loc[bg.drug_connected & in_seed, "ensembl_id"]) == set(edges.Gene_ID) & set(s_in.ensembl_id), \
        "any-drug flags disagree with Task 2"
    assert set(bg.loc[bg.pharm_connected & in_seed, "ensembl_id"]) == \
        set(edges.loc[edges.pharmacological, "Gene_ID"]) & set(s_in.ensembl_id), "pharm flags disagree with Task 2"

    groups = [("All seed genes", "summary", s_in.ensembl_id),
              ("Histology genes (stages 1-4)", "summary", s_in.loc[s_in.main_stage != "0_genetic_risk", "ensembl_id"])]
    for st in STAGES:
        groups.append((st, "stage", s_in.loc[s_in.main_stage == st, "ensembl_id"]))
    for st in STAGES[1:]:
        for d in ["increasing", "decreasing"]:
            ids = s_in.loc[(s_in.main_stage == st) & (s_in.direction == d), "ensembl_id"]
            if len(ids):
                groups.append((f"{st}_{d}", "stage x direction", ids))

    bins = bg.degree_bin.values
    bin_index = {b: np.flatnonzero(bins == b) for b in np.unique(bins)}

    def analyse(col, et):
        rng = np.random.default_rng(SEED)
        conn = bg[col].values
        dpg = et.groupby("Gene_ID").Drug_ID.agg(set)
        rows = []
        for name, level, ids in groups:
            ids = pd.Index(ids)
            sub = bgi.loc[ids]
            k, n = int(sub[col].sum()), len(sub)
            rest = bgi.drop(ids)
            k0, n0 = int(rest[col].sum()), len(rest)
            table = [[k, n - k], [k0, n0 - k0]]
            orr = odds_ratio(table)
            ci = orr.confidence_interval(0.95)
            lo, hi = wilson(k, n)
            need = sub.degree_bin.value_counts()
            null = np.empty(N_PERM)
            for i in range(N_PERM):
                pick = np.concatenate([rng.choice(bin_index[b], size=c, replace=False) for b, c in need.items()])
                null[i] = conn[pick].mean()
            obs = k / n
            p_enr = (1 + (null >= obs).sum()) / (N_PERM + 1)
            p_dep = (1 + (null <= obs).sum()) / (N_PERM + 1)
            rows.append({"Group": name, "Level": level, "Total genes": n, "Drug-connected genes": k,
                         "Percentage": round(100 * obs, 2), "Wilson_CI_low_%": round(100 * lo, 2),
                         "Wilson_CI_high_%": round(100 * hi, 2),
                         "Unique drugs": len(set().union(*[dpg.get(g, set()) for g in ids])),
                         "Median PPI degree": float(sub.ppi_degree.median()),
                         "Rest-of-background %": round(100 * k0 / n0, 2), "Odds ratio": round(orr.statistic, 3),
                         "OR_CI_low": round(ci.low, 3), "OR_CI_high": round(ci.high, 3),
                         "Fisher_p": fisher_exact(table).pvalue,
                         "Degree-matched expected %": round(100 * null.mean(), 2),
                         "Degree-matched null 95% range %": f"{100 * np.quantile(null, 0.025):.1f}-{100 * np.quantile(null, 0.975):.1f}",
                         "Fold vs degree-matched": round(obs / null.mean(), 3) if null.mean() else np.nan,
                         "Perm_p_enrichment": p_enr, "Perm_p_depletion": p_dep,
                         "Perm_p_two_sided": min(1.0, 2 * min(p_enr, p_dep))})
        df = pd.DataFrame(rows)
        df["Fisher_BH_FDR"] = false_discovery_control(df.Fisher_p, method="bh")
        df["Perm_BH_FDR"] = false_discovery_control(df.Perm_p_two_sided, method="bh")
        pairs = []
        for a, b in combinations(STAGES[1:], 2):
            A, B = bgi.loc[s_in.loc[s_in.main_stage == a, "ensembl_id"]], bgi.loc[s_in.loc[s_in.main_stage == b, "ensembl_id"]]
            ka, na, kb, nb = int(A[col].sum()), len(A), int(B[col].sum()), len(B)
            t = [[ka, na - ka], [kb, nb - kb]]
            r = odds_ratio(t)
            ci = r.confidence_interval(0.95)
            pairs.append({"Comparison": f"{a} vs {b}", "%A": round(100 * ka / na, 1), "%B": round(100 * kb / nb, 1),
                          "Odds ratio": round(r.statistic, 3), "OR_CI_low": round(ci.low, 3),
                          "OR_CI_high": round(ci.high, 3), "Fisher_p": fisher_exact(t).pvalue})
        pairs = pd.DataFrame(pairs)
        pairs["BH_FDR"] = false_discovery_control(pairs.Fisher_p, method="bh")
        return df, pairs

    res = {}
    for key, (col, label, suffix) in DEFS.items():
        df, pairs = analyse(col, edges if key == "any" else edges[edges.pharmacological])
        df.to_csv(OUT / f"druggability_summary{suffix}.csv", index=False)
        pairs.to_csv(OUT / f"pairwise_stage_comparisons{suffix}.csv", index=False)
        res[key] = (df, pairs)

    fig, axes = plt.subplots(2, 1, figsize=(13, 10))
    for ax, key in zip(axes, ["pharm", "any"]):
        col, label, _ = DEFS[key]
        df = res[key][0]
        x = np.arange(len(df))
        for xi, (_, r) in zip(x, df.iterrows()):
            a = 1.0 if r["Level"] != "stage x direction" else 0.55
            ax.bar(xi, r["Percentage"], color="#2b6cb0", alpha=a)
            ax.errorbar(xi, r["Percentage"], yerr=[[r["Percentage"] - r["Wilson_CI_low_%"]],
                                                    [r["Wilson_CI_high_%"] - r["Percentage"]]], color="#1f2933", capsize=3, lw=1)
            ax.plot([xi - 0.4, xi + 0.4], [r["Degree-matched expected %"]] * 2, color="#1f2933", ls=":", lw=1.4)
            ax.text(xi, r["Wilson_CI_high_%"] + 1, f"{r['Percentage']:.0f}%{'*' if r['Perm_BH_FDR'] < 0.05 else ''}",
                    ha="center", fontsize=8)
        bgp = 100 * bg[col].mean()
        ax.axhline(bgp, color="#616e7c", ls="--", lw=1.2)
        ax.text(len(df) - 0.5, bgp + 0.8, f"Background {bgp:.1f}%", ha="right", fontsize=9, color="#616e7c")
        ax.set_xticks(x, [f"{g}\n(n={n})" for g, n in zip(df.Group, df["Total genes"])], fontsize=7.5, rotation=30, ha="right")
        ax.set_ylabel("% genes drug-connected")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"{'A' if key == 'pharm' else 'B'}. {label}", loc="left", fontweight="bold")
    fig.text(0.07, 0.005, "Bars: observed % (95% Wilson CI). Dotted: degree-matched expected % (1,000 permutations, seed 42). "
                          "* = permutation BH-FDR < 0.05.", fontsize=8.5, color="#616e7c")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(OUT / "druggability_barplot.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    deg = bg.groupby("degree_bin").agg(genes=("ensembl_id", "size"), min_degree=("ppi_degree", "min"),
                                       max_degree=("ppi_degree", "max"), pct_any=("drug_connected", "mean"),
                                       pct_pharm=("pharm_connected", "mean"))
    deg[["pct_any", "pct_pharm"]] = (100 * deg[["pct_any", "pct_pharm"]]).round(1)
    deg.to_csv(OUT / "background_degree_bins.csv")

    def md(t, cols):
        t = t[cols]
        return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
                         + ["| " + " | ".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in r) + " |" for r in t.values])

    cols = ["Group", "Total genes", "Drug-connected genes", "Percentage", "Degree-matched expected %",
            "Fold vs degree-matched", "Odds ratio", "Fisher_BH_FDR", "Perm_BH_FDR"]
    summary = f"""# Task 3: druggability of the seed genes

Background: **{len(bg):,} genes** with at least one INTERACTS_WITH partner in this graph
(pharmacological target {100 * bg.pharm_connected.mean():.1f}%, any drug {100 * bg.drug_connected.mean():.1f}%).
{len(s_in)} of 280 seed genes are in the background; {len(outside)} have no partner and are left out
(`seed_genes_outside_background.csv`). Drug flags for seed genes match Task 2. Generated {now()}.

## A. Pharmacological target (main definition)

{md(res['pharm'][0], cols)}

Pairwise stage comparisons (Fisher, BH-FDR): see `pairwise_stage_comparisons_pharmacological.csv`.

{md(res['pharm'][1], ['Comparison', '%A', '%B', 'Odds ratio', 'Fisher_p', 'BH_FDR'])}

## B. Any drug link (secondary)

{md(res['any'][0], cols)}

Definitions: pharmacological target = an OptimusKG action other than ENZYME / TRANSPORTER / CARRIER, a curated link added in step 0.4,
or an OptimusKG link carrying a step 0.5 direction. Fisher's exact test compares a group with the rest of the background.
The permutation test draws 1,000 random background gene sets matched on INTERACTS_WITH degree decile (seed 42); well-studied genes
have more partners and more drugs, so the degree-matched comparison is the fairer one.

Outputs: `druggability_summary_pharmacological.csv`, `druggability_summary_any.csv`, `pairwise_stage_comparisons_*.csv`,
`background_genes_druggability.csv`, `background_degree_bins.csv`, `druggability_barplot.png`. Script: `code/task3_druggability.py`.
"""
    (OUT / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
