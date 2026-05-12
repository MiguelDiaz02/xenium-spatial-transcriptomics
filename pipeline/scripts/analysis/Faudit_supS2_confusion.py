#!/usr/bin/env python3
"""
Faudit_supS2_confusion.py — Supplementary Figure S2 (Banksy × Novae confusion matrix).

Renders the cell-level confusion matrix between the two spatial domain partitionings
referenced in §3.3 of the audited manuscript. Shows where the 35.4% consensus comes
from and where the 64.6% boundary/transition cells live.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
DOM_DIR = BASE / "human_lung_cancer/results/03_phase3_spatial/F1_spatial_domains"
OUT_MS = BASE / "manuscript/figures/FigS2_banksy_novae_confusion.png"
OUT_RES = BASE / "human_lung_cancer/results/figures/ms_sup_figS2_banksy_novae_confusion.png"

print("Loading per-cell domain assignments ...")
banksy = pd.read_csv(DOM_DIR / "banksy_domains.csv")
novae = pd.read_csv(DOM_DIR / "novae_domains.csv")
print(f"  Banksy:  shape={banksy.shape}, cols={list(banksy.columns)[:5]}")
print(f"  Novae:   shape={novae.shape}, cols={list(novae.columns)[:5]}")

# Use first column as cell ID, second as label
b_id_col = banksy.columns[0]
n_id_col = novae.columns[0]
b_lab_col = "banksy_domain" if "banksy_domain" in banksy.columns else \
            [c for c in banksy.columns if "domain" in c.lower() or "cluster" in c.lower()][0]
# Prefer the canonical novae_domain column (not novae_domains_LX multilevel)
if "novae_domain" in novae.columns:
    n_lab_col = "novae_domain"
else:
    n_lab_col = [c for c in novae.columns if "domain" in c.lower()][0]
print(f"  Banksy label col: '{b_lab_col}'   Novae label col: '{n_lab_col}'")

merged = banksy[[b_id_col, b_lab_col]].merge(
    novae[[n_id_col, n_lab_col]], left_on=b_id_col, right_on=n_id_col, how="inner"
)
merged = merged[[b_lab_col, n_lab_col]].dropna()
merged.columns = ["banksy", "novae"]
print(f"  Merged cells with both labels: {len(merged):,}")

ct = pd.crosstab(merged["banksy"], merged["novae"], normalize="index") * 100
ct_counts = pd.crosstab(merged["banksy"], merged["novae"])

ari = adjusted_rand_score(merged["banksy"].astype(str), merged["novae"].astype(str))
nmi = normalized_mutual_info_score(merged["banksy"].astype(str), merged["novae"].astype(str))
print(f"  ARI(Banksy, Novae)  = {ari:.4f}")
print(f"  NMI(Banksy, Novae)  = {nmi:.4f}")

# Sort Banksy rows by best-matched Novae column (preserves block structure;
# avoids the duplicate-column artifact of a greedy reuse strategy).
best_novae = ct.idxmax(axis=1)
banksy_order = best_novae.sort_values().index.tolist()
ct = ct.loc[banksy_order]
ct_counts = ct_counts.loc[banksy_order]
# Sort Novae columns by total cell count (descending)
col_totals = ct_counts.sum(axis=0).sort_values(ascending=False)
ct = ct.reindex(columns=col_totals.index)
ct_counts = ct_counts.reindex(columns=col_totals.index)

fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="white")

# Panel A — confusion (% per Banksy row)
ax = axes[0]
im = ax.imshow(ct.values, cmap="YlOrRd", aspect="auto", vmin=0, vmax=100)
ax.set_xticks(range(len(ct.columns)))
ax.set_xticklabels([str(c) for c in ct.columns], rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(ct.index)))
ax.set_yticklabels([str(c) for c in ct.index], fontsize=9)
ax.set_xlabel("Novae domain")
ax.set_ylabel("Banksy domain")
ax.set_title(f"A  Cell assignment overlap (% per Banksy domain)\n"
             f"ARI = {ari:.3f}   NMI = {nmi:.3f}",
             loc="left", fontsize=11, fontweight="bold")
plt.colorbar(im, ax=ax, label="% of Banksy domain", shrink=0.8)

# Annotate cells with values >= 5%
for i in range(ct.shape[0]):
    for j in range(ct.shape[1]):
        v = ct.values[i, j]
        if v >= 5:
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    fontsize=7, color="white" if v > 50 else "black")

# Panel B — alignment quality histogram
ax2 = axes[1]
diag_pct = []
for i, r in enumerate(ct.index):
    if i < ct.shape[1]:
        diag_pct.append(ct.iloc[i, i])
ax2.hist(diag_pct, bins=10, color="#4575b4", alpha=0.85, edgecolor="white")
ax2.axvline(np.mean(diag_pct), color="#d73027", linestyle="--", lw=1.5,
            label=f"mean = {np.mean(diag_pct):.1f}%")
ax2.set_xlabel("Best-aligned overlap (%)", fontsize=11)
ax2.set_ylabel("Number of Banksy domains", fontsize=11)
ax2.set_title("B  Per-Banksy-domain best alignment with Novae\n"
              "(higher = more agreement on that domain's cell membership)",
              loc="left", fontsize=11, fontweight="bold")
ax2.legend()

plt.tight_layout()
for outpath in (OUT_MS, OUT_RES):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"Saved: {outpath}")
plt.close()

# Save summary
summary_path = BASE / "human_lung_cancer/results/02_biology/audit_corrected/banksy_novae_confusion_summary.csv"
summary_path.parent.mkdir(parents=True, exist_ok=True)
ct_counts.to_csv(summary_path)
print(f"\nConfusion (counts): {summary_path}")
print(f"ARI = {ari:.4f}, NMI = {nmi:.4f}")
