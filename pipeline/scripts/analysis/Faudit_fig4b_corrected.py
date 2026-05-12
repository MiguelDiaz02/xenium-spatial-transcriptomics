#!/usr/bin/env python3
"""
Faudit_fig4b_corrected.py — regenerate Fig4B with Bonferroni-corrected interactions.

Replaces the OLD Fig4B which showed 22 dual-validated interactions (per-test).
The audited count is 20 of 30 submitted (Bonferroni at p < 1.67e-3).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
SPACIA_RES = BASE / "human_lung_cancer/results/02_biology/audit_corrected/spacia_results_corrected.csv"
OUT_MS = BASE / "manuscript/figures/Fig4B_spacia_axes.png"
OUT_RES = BASE / "human_lung_cancer/results/figures/ms_main_fig4B_spacia_axes.png"

df = pd.read_csv(SPACIA_RES)
df = df[df["pathway_pval"].notna()].copy()
n_total_attempted = 30  # 5 failed + 25 completed
bonf_thresh = 0.05 / n_total_attempted

df["axis"] = df["ligand_complex"] + " → " + df["receptor_complex"]
df["validated"] = df["pathway_pval"] < bonf_thresh
df["neg_log_p"] = -np.log10(df["pathway_pval"].clip(1e-300))

# Color per axis
axis_colors = {
    "ADAM17 → MUC1": "#d62728",
    "LTF → AGER": "#1f77b4",
    "CDH1 → EGFR": "#2ca02c",
    "S100B → AGER": "#9467bd",
    "CD24 → SELP": "#bcbd22",
}
df["color"] = df["axis"].map(lambda a: axis_colors.get(a, "#888888"))

fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), facecolor="white")

# ── Panel A — scatter ─────────────────────────────────────────────────────
ax = axes[0]
val_mask = df["validated"]
ax.scatter(df.loc[val_mask, "scaled_weight"], df.loc[val_mask, "neg_log_p"],
           c=df.loc[val_mask, "color"], s=80, alpha=0.85,
           edgecolors="white", linewidths=0.6, zorder=3, label="_nolegend_")
ax.scatter(df.loc[~val_mask, "scaled_weight"], df.loc[~val_mask, "neg_log_p"],
           c="#cccccc", s=60, alpha=0.7, marker="x", zorder=2, label="_nolegend_")
ax.axhline(-np.log10(bonf_thresh), color="#444", linestyle="--", lw=1.2,
           label=f"Bonferroni (p<{bonf_thresh:.2e})")
ax.axhline(-np.log10(0.05), color="#888", linestyle=":", lw=1.0,
           label="Per-test (p<0.05)")
ax.set_xlabel("LIANA+ scaled interaction score", fontsize=11)
ax.set_ylabel(r"$-\log_{10}$($p_{\text{Spacia}}$)", fontsize=11)
ax.set_title(f"A  Bonferroni-corrected dual-validated interactions\n"
             f"({val_mask.sum()} of {n_total_attempted} submitted; "
             f"{(~val_mask).sum()} per-test only)",
             loc="left", fontsize=11, fontweight="bold")
legend_axes = [Patch(color=c, label=a) for a, c in axis_colors.items()
               if a in df["axis"].values]
legend_axes.append(Patch(color="#cccccc", label="Not Bonf-validated"))
ax.legend(handles=legend_axes, loc="upper left", fontsize=8, ncol=1)

# ── Panel B — bars per axis with sender count ─────────────────────────────
ax2 = axes[1]
axis_summary = (df[df["validated"]]
                .groupby("axis")
                .agg(n_unique_senders=("source", "nunique"),
                     n_pairs=("source", "size"),
                     max_neg_log_p=("neg_log_p", "max"))
                .sort_values("max_neg_log_p", ascending=True))

bar_colors = [axis_colors.get(a, "#888888") for a in axis_summary.index]
y = np.arange(len(axis_summary))
ax2.barh(y, axis_summary["max_neg_log_p"], color=bar_colors, alpha=0.85,
         edgecolor="white", linewidth=1)
ax2.set_yticks(y)
ax2.set_yticklabels(axis_summary.index, fontsize=10)
for i, (axis, row) in enumerate(axis_summary.iterrows()):
    ax2.text(row["max_neg_log_p"] + 1.5, i,
             f"{int(row['n_unique_senders'])} senders / "
             f"{int(row['n_pairs'])} pairs",
             va="center", fontsize=9, color="#222")
ax2.set_xlabel(r"max $-\log_{10}$($p_{\text{Spacia}}$) per axis", fontsize=11)
ax2.set_title("B  Bonferroni-validated signaling axes\n"
              "(unique senders × validated sender→receiver pairs)",
              loc="left", fontsize=11, fontweight="bold")
ax2.set_xlim(0, axis_summary["max_neg_log_p"].max() * 1.25)

plt.tight_layout()
for outpath in (OUT_MS, OUT_RES):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"Saved: {outpath}")
plt.close()

print("\nValidated axes summary:")
print(axis_summary.to_string())
