"""Figures of the evaluation (short labels, large text): robustness curves and evidence composition.

    python code/eval_figures.py [output_dir]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "results/Evaluation"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else EV / "figures"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 15, "axes.labelsize": 16, "axes.titlesize": 16, "legend.fontsize": 13,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.spines.top": False, "axes.spines.right": False})
METHODS = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]
COL = {"SP": "#9aa0a6", "NP": "#7aa6c2", "RWR": "#4f81a8", "MPC": "#e0a458", "DWPC": "#b9822c", "HeteSim": "#8a6bb1", "Proposed": "#c0392b"}
MARK = {"SP": "x", "NP": "s", "RWR": "^", "MPC": "o", "DWPC": "D", "HeteSim": "v", "Proposed": "o"}


def robustness():
    r = pd.read_csv(EV / "robust/robustness_long.csv")
    r = r[(r.K == 10) & r.kind.isin(["seeds", "edges", "noop"])].copy()
    r["level"] = r.level.astype(float)
    ceil = r[r.kind == "noop"].groupby("method").family_jaccard.mean()
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
    for a, kind, title in zip(ax, ["seeds", "edges"], ["Seeds removed", "Edges removed"]):
        d = r[r.kind == kind].groupby(["level", "method"]).family_jaccard.mean().unstack()
        x = [0] + list(d.index * 100)
        for m in METHODS:
            if m == "SP":
                continue
            y = [ceil[m]] + list(d[m])
            a.plot(x, y, marker=MARK[m], color=COL[m], lw=3.2 if m == "Proposed" else 2, ms=7, label=m)
        a.set_title(title)
        a.set_xlabel("% removed")
        a.set_xticks(x)
        a.set_ylim(0, 1.03)
        a.grid(alpha=0.25)
    ax[0].set_ylabel("Top-10 overlap")
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, frameon=False, columnspacing=1.2, handlelength=1.6)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUT / "g_robust.png", dpi=220)
    plt.close(fig)


def evidence_bars():
    t = pd.read_csv(EV / "tables/evidence_all_configs.csv")
    t = t[(t.K == 20) & t.method.isin(METHODS) & (t.method != "SP")].set_index("method").reindex([m for m in METHODS if m != "SP"])
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    bottom = np.zeros(len(t))
    cols = {"human": "#2e7d32", "animal": "#81b29a", "in_vitro": "#cfe3d4"}
    for k, lab in [("human", "Human"), ("animal", "Animal"), ("in_vitro", "Cells")]:
        ax.bar(t.index, t[k], bottom=bottom, color=cols[k], label=lab, edgecolor="white")
        bottom += t[k].to_numpy()
    ax.bar(t.index, 20 - bottom, bottom=bottom, color="#e6e6e6", label="None", edgecolor="white")
    ax.set_ylabel("Drugs in top 20")
    ax.set_ylim(0, 20)
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.2), frameon=False, columnspacing=1.0, handlelength=1.2)
    plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
    fig.tight_layout()
    fig.savefig(OUT / "g_evidence.png", dpi=220)
    plt.close(fig)


def heatmap():
    o = pd.read_csv(EV / "main/overlap_jaccard.csv")
    o = o[o.K == 20].pivot(index="a", columns="b", values="jaccard").reindex(index=METHODS, columns=METHODS)
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    im = ax.imshow(o.to_numpy(), cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(METHODS)), METHODS, rotation=35, ha="right")
    ax.set_yticks(range(len(METHODS)), METHODS)
    for i in range(len(METHODS)):
        for j in range(len(METHODS)):
            v = o.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=11, color="white" if v > 0.55 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.tight_layout()
    fig.savefig(OUT / "g_overlap_methods.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    heatmap()
    robustness()
    if (EV / "tables/evidence_all_configs.csv").exists():
        evidence_bars()
    print("figures in", OUT)
