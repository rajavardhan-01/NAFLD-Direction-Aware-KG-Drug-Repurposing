"""Stability measured from each method's own unperturbed value.

Tied scores are broken at random, so a ranking compared with itself (no perturbation, "noop") can overlap by less than 1.
At target-family level this value is below 1 for the proposed method only, among the methods compared here. For each
perturbation level this script computes, for every method, the loss from its own unperturbed value, and the paired
difference between the loss of the proposed method and the loss of each baseline over the same perturbed graphs
(mean with a normal-approximation 95% interval, mean +/- 1.96 standard errors). The unperturbed value of a method is
its mean over the noop repetitions and is treated as a constant, so the intervals do not include its sampling error.

SP is left out: 233 of 299 drugs share its top score.

Reads results/Evaluation/robust/robustness_long.csv (written by eval_robust.py); no database is needed.

    python code/eval_tie_adjusted.py     # writes results/Evaluation/robust/tie_adjusted_stability.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ROB = ROOT / "results/Evaluation/robust"
OUT = ROB / "tie_adjusted_stability.csv"

METHODS = ["NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]
K = 10


def main():
    long = pd.read_csv(ROB / "robustness_long.csv")
    long = long[(long.K == K) & long.method.isin(METHODS)].copy()
    long["level"] = long.level.astype(str)
    noop = long[long.kind == "noop"].groupby("method")[["family_jaccard", "jaccard"]].mean()
    rows = []
    runs = long[long.kind.isin(["seeds", "edges"])]
    for (kind, level), group in runs.groupby(["kind", "level"]):
        overlap = group.pivot_table(index="rep", columns="method", values="family_jaccard")
        loss = noop.family_jaccard[METHODS] - overlap[METHODS]
        for m in METHODS:
            row = {"kind": kind, "level": float(level), "method": m, "reps": int(loss[m].count()),
                   "unperturbed_drug_overlap": noop.jaccard[m], "unperturbed_family_overlap": noop.family_jaccard[m],
                   "family_overlap": overlap[m].mean(), "loss": loss[m].mean(),
                   "proposed_minus_method_loss": np.nan, "ci_low": np.nan, "ci_high": np.nan, "reading": ""}
            if m != "Proposed":
                diff = (loss["Proposed"] - loss[m]).dropna()
                half = 1.96 * diff.std() / np.sqrt(len(diff))
                low, high = diff.mean() - half, diff.mean() + half
                row.update(proposed_minus_method_loss=diff.mean(), ci_low=low, ci_high=high,
                           reading="proposed loses more" if low > 0 else "proposed loses less" if high < 0 else "no clear difference")
            rows.append(row)
    out = pd.DataFrame(rows).sort_values(["kind", "level", "method"], ascending=[False, True, True])
    out.round(4).to_csv(OUT, index=False)
    print("unperturbed family overlap:", noop.family_jaccard[METHODS].round(3).to_dict())
    print(out[out.method != "Proposed"].pivot_table(index="method", columns="reading", values="reps", aggfunc="size", fill_value=0).to_string())
    print("written", OUT)


if __name__ == "__main__":
    main()
