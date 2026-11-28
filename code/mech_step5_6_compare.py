"""Steps 5 and 6: compare retained / excluded / conventional-MPC / annotation-only drugs on mechanistic evidence (M0-M4),
with confidence intervals, contradictory evidence and target-family sensitivity analyses.

Drug-level grade = best M over the drug's seed-gene pairs (the pairs that start its routes). Groups (they overlap):
  retained_38          the 38 drugs ranked by the direction-aware method
  excluded_78          the 78 further drugs of the conventional MPC list
  mpc_116              all 116 conventional-MPC drugs (= retained + excluded)
  annotation_only_39   drugs ranked when the filter only needs a known direction (decidability-only filter)
Target family = the set of seed targets of a drug; drugs with the same set act through one mechanism.

Outputs in results/Evaluation/mechanistic/: step5_drug_level.csv, step5_group_comparison.csv, step6_difference_tests.csv,
step6_family_sensitivity.csv, step6_contradictions.csv and summary.md

    python code/mech_step5_6_compare.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/Evaluation/mechanistic"
rng = np.random.default_rng(20261122)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def boot_mean(x, B=10000):
    x = np.asarray(x, float)
    m = rng.choice(x, size=(B, len(x))).mean(1)
    return np.percentile(m, [2.5, 97.5])


groups = pd.read_csv(OUT / "step1_drug_groups.csv")
pg = pd.read_csv(OUT / "step4_pair_grades.csv")
# contradictory evidence: automatic flags were re-read with an AI assistant, not verified by a human expert
# (decisions in step4_manual_review.csv; the file name is historical); only confirmed ones count
rev = pd.read_csv(OUT / "step4_manual_review.csv")
conf = set(zip(rev[rev.decision == "confirmed"].Drug, rev[rev.decision == "confirmed"].Gene))
pg["contradiction_auto"] = pg.contradictory_evidence
pg["contradictory_evidence"] = [(a, b) in conf for a, b in zip(pg.Drug, pg.Gene)]
pg.to_csv(OUT / "step4_pair_grades.csv", index=False)
d = pg.groupby("Drug_ID").agg(auto_flag=("contradiction_auto", "max"), M=("M", "max"), pairs=("Gene", "nunique"), contradictory=("contradictory_evidence", "max"),
                              opposite=("opposite_only", "sum"), negative=("negative", "sum"),
                              targets=("Gene", lambda s: ";".join(sorted(set(s))))).reset_index()
d = groups[["Drug_ID", "Drug", "group", "annotation_only_group", "retained_by_filter", "mpc_rank"]].merge(d, on="Drug_ID", how="left")
d["family"] = d.targets
d.to_csv(OUT / "step5_drug_level.csv", index=False)
masks = {"retained_38": d.retained_by_filter, "excluded_78": ~d.retained_by_filter, "mpc_116": d.retained_by_filter | ~d.retained_by_filter,
         "annotation_only_39": d.annotation_only_group}


def summarise(x, unit="drug"):
    n = len(x)
    row = {"unit": unit, "n": n}
    for lab, thr in (("M>=2", 2), ("M>=3", 3), ("M4", 4)):
        k = int((x.M >= thr).sum()); lo, hi = wilson(k, n)
        row[lab] = f"{k}/{n} ({100 * k / n:.0f}%)"; row[lab + "_ci"] = f"{100 * lo:.0f}-{100 * hi:.0f}%"
    lo, hi = boot_mean(x.M)
    row["mean_M"] = round(x.M.mean(), 2); row["mean_M_ci"] = f"{lo:.2f}-{hi:.2f}"
    k = int(x.contradictory.sum()); lo, hi = wilson(k, n)
    row["contradictory"] = f"{k}/{n} ({100 * k / n:.0f}%)"; row["contradictory_ci"] = f"{100 * lo:.0f}-{100 * hi:.0f}%"
    if "auto_flag" in x:
        row["auto_flagged_before_review"] = int(x.auto_flag.sum())
    row["M0"] = int((x.M == 0).sum()); row["M1"] = int((x.M == 1).sum()); row["M2"] = int((x.M == 2).sum()); row["M3"] = int((x.M == 3).sum())
    return row


comp = pd.DataFrame([{"group": g, **summarise(d[m])} for g, m in masks.items()])
comp.to_csv(OUT / "step5_group_comparison.csv", index=False)

# ---- difference tests (drug level)
tests = []
for a, b in (("retained_38", "excluded_78"), ("annotation_only_39", "retained_38"), ("retained_38", "mpc_116")):
    A_, B_ = d[masks[a]], d[masks[b]]
    for lab, thr in (("M>=2", 2), ("M>=3", 3), ("M4", 4)):
        ka, kb = int((A_.M >= thr).sum()), int((B_.M >= thr).sum())
        diff = ka / len(A_) - kb / len(B_)
        bs = [rng.choice((A_.M >= thr).values, len(A_)).mean() - rng.choice((B_.M >= thr).values, len(B_)).mean() for _ in range(10000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        p = fisher_exact([[ka, len(A_) - ka], [kb, len(B_) - kb]])[1] if a != "retained_38" or b != "mpc_116" else np.nan
        tests.append({"comparison": f"{a} vs {b}", "threshold": lab, "share_a": round(ka / len(A_), 3), "share_b": round(kb / len(B_), 3),
                      "difference": round(diff, 3), "boot_95ci": f"{lo:+.2f} to {hi:+.2f}", "fisher_p": None if np.isnan(p) else round(p, 3)})
pd.DataFrame(tests).to_csv(OUT / "step6_difference_tests.csv", index=False)

# ---- target-family sensitivity
fam_rows = []
for g, m in masks.items():
    x = d[m]
    fam = x.groupby("family").agg(M=("M", "max"), contradictory=("contradictory", "max"), drugs=("Drug_ID", "nunique")).reset_index()
    fam_rows.append({"group": g, "analysis": "all drugs (drug unit)", **summarise(x)})
    fam_rows.append({"group": g, "analysis": "one value per target family (best M)", **summarise(fam, "family")})
    top = fam.sort_values("drugs", ascending=False).iloc[0].family
    fam_rows.append({"group": g, "analysis": f"drug unit without the largest family ({top})", **summarise(x[x.family != top])})
    fam_rows.append({"group": g, "analysis": "drug unit without PPARG/PPARA/TNF families", **summarise(x[~x.family.str.contains(r"\bPPARG\b|\bPPARA\b|\bTNF\b")])})
pd.DataFrame(fam_rows).to_csv(OUT / "step6_family_sensitivity.csv", index=False)

# ---- contradictory evidence list
c = pg[pg.contradictory_evidence].merge(groups[["Drug_ID", "group"]], on="Drug_ID")
c[["group", "Drug", "Gene", "annotated_direction", "M", "opposite_only", "negative", "pmids_contradiction"]].to_csv(OUT / "step6_contradictions.csv", index=False)

pd.set_option("display.width", 250)
print(comp[["group", "n", "M>=2", "M>=3", "M4", "mean_M", "contradictory"]].to_string(index=False))
print(pd.DataFrame(tests).to_string(index=False))
print(pd.DataFrame(fam_rows)[["group", "analysis", "n", "M>=3", "M>=3_ci", "M4"]].to_string(index=False))
