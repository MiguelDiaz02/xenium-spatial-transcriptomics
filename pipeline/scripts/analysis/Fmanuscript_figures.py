"""
Generate all manuscript figures not yet produced.
Figures: 2A (UMAP), 2B (spatial cell types), 2C (ARI summary),
         3A/3B/3C (spatial domains), 4B (Spacia axes), 5B (pseudotime proportions),
         S1 (QC metrics).
Environment: xenium_pipeline
"""

import os
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgba
import matplotlib.gridspec as gridspec

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import dataset_root

BASE = str(dataset_root())
RESULTS = os.path.join(BASE, "results")
SPATIAL = os.path.join(RESULTS, "03_phase3_spatial")
FIG_OUT = os.path.join(RESULTS, "figures", "manuscript")
os.makedirs(FIG_OUT, exist_ok=True)

DPI = 300
PANEL_FONT = 14

# ── Crameri / fallback palettes ───────────────────────────────────────────────
try:
    from cmcrameri import cm as cmc
    CMAP_SEQ  = cmc.batlow
    CMAP_DIV  = cmc.vik
    print("Using Crameri colormaps.")
except ImportError:
    CMAP_SEQ  = plt.cm.viridis
    CMAP_DIV  = plt.cm.RdBu_r
    print("cmcrameri not found; using matplotlib defaults.")

# 26 distinct colors for cell types (colorblind-aware)
CELL_PALETTE = [
    "#4C72B0","#DD8452","#55A868","#C44E52","#8172B2",
    "#937860","#DA8BC3","#8C8C8C","#CCB974","#64B5CD",
    "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
    "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF",
    "#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5",
    "#C49C94"
]

def label_panel(ax, letter, fontsize=16):
    ax.text(-0.08, 1.05, letter, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="right")

def save(fig, name):
    path = os.path.join(FIG_OUT, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {name}")


# ═══════════════════════════════════════════════════════════════════════════════
# LOAD SDATA (obs only — avoids loading expression matrix)
# ═══════════════════════════════════════════════════════════════════════════════
print("Loading sdata (obs + obsm only)...")
import spatialdata as sd
sdata = sd.read_zarr(os.path.join(RESULTS, "sdata.zarr"))
adata = sdata.tables["table"]
obs   = adata.obs.copy()
xy    = adata.obsm["spatial"]          # (268034, 2)  x, y
umap  = adata.obsm["X_umap"]           # (268034, 2)
cell_types_ordered = obs["cell_type_L2"].value_counts().index.tolist()
ct_color = {ct: CELL_PALETTE[i % len(CELL_PALETTE)]
            for i, ct in enumerate(cell_types_ordered)}
print(f"  {len(obs)} cells, {obs['cell_type_L2'].nunique()} L2 types.")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 2A — UMAP colored by cell_type_L2 with abundance labels
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 2A — UMAP...")
fig, ax = plt.subplots(figsize=(9, 7))

pct = obs["cell_type_L2"].value_counts(normalize=True) * 100
# plot each type
for ct in cell_types_ordered:
    mask = obs["cell_type_L2"] == ct
    ax.scatter(umap[mask, 0], umap[mask, 1],
               c=ct_color[ct], s=0.3, alpha=0.5, linewidths=0, rasterized=True)

handles = [mpatches.Patch(color=ct_color[ct],
                           label=f"{ct}  ({pct.get(ct,0):.1f}%)")
           for ct in cell_types_ordered]
ax.legend(handles=handles, fontsize=6, loc="upper right",
          framealpha=0.9, ncol=2, markerscale=3,
          handlelength=1.2, borderpad=0.8)
ax.set_xlabel("UMAP 1", fontsize=11)
ax.set_ylabel("UMAP 2", fontsize=11)
ax.set_title("Cell-type annotation — 26 populations", fontsize=13)
ax.axis("off")
label_panel(ax, "A")
fig.tight_layout()
save(fig, "Fig2A_umap_celltypes.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 2B — Spatial XY map colored by cell_type_L2
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 2B — Spatial cell-type map...")
fig, ax = plt.subplots(figsize=(9, 8))

for ct in cell_types_ordered:
    mask = (obs["cell_type_L2"] == ct).values
    ax.scatter(xy[mask, 0], xy[mask, 1],
               c=ct_color[ct], s=0.15, alpha=0.6, linewidths=0, rasterized=True)

handles = [mpatches.Patch(color=ct_color[ct], label=ct)
           for ct in cell_types_ordered]
ax.legend(handles=handles, fontsize=6, loc="upper right",
          framealpha=0.9, ncol=2, markerscale=3,
          handlelength=1.2, borderpad=0.8)
ax.set_xlabel("X coordinate (µm)", fontsize=11)
ax.set_ylabel("Y coordinate (µm)", fontsize=11)
ax.set_title("Spatial distribution of 26 cell populations", fontsize=13)
ax.set_aspect("equal")
label_panel(ax, "B")
fig.tight_layout()
save(fig, "Fig2B_spatial_celltypes.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 2C — Annotation validation: ARI bar chart + source breakdown
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 2C — Annotation ARI summary...")
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Panel left: ARI values per comparison
comparisons = [
    "vs. immune\nsubclustering",
    "vs. Leiden\nclusters",
    "vs. Leiden\n(immune only)",
]
aris = [0.742, 0.293, 0.194]
colors_bar = ["#2196F3", "#90CAF9", "#BBDEFB"]
bars = axes[0].barh(comparisons, aris, color=colors_bar, edgecolor="white", height=0.55)
for bar, val in zip(bars, aris):
    axes[0].text(val + 0.01, bar.get_y() + bar.get_height()/2,
                 f"{val:.3f}", va="center", fontsize=11, fontweight="bold")
axes[0].set_xlim(0, 0.90)
axes[0].set_xlabel("Adjusted Rand Index (ARI)", fontsize=11)
axes[0].set_title("Annotation agreement\nwith orthogonal references", fontsize=11)
axes[0].axvline(0.742, color="#1565C0", lw=1.2, ls="--", alpha=0.5)
axes[0].spines[["top","right"]].set_visible(False)
label_panel(axes[0], "C")

# Panel right: annotation source breakdown
src_counts = obs["annotation_source_v3"].value_counts()
src_labels = {
    "immune_granular":    "Immune\nsubclustering",
    "score_nonimmune":    "Marker\nscoring",
}
# handle dynamic leiden keys
for k in src_counts.index:
    if k.startswith("leiden") and k not in src_labels:
        src_labels[k] = "Leiden\ndirect"

labels_plot = [src_labels.get(k, k) for k in src_counts.index]
wedge_colors = ["#4CAF50", "#FF9800", "#9C27B0", "#F44336"][:len(src_counts)]
wedges, texts, autotexts = axes[1].pie(
    src_counts.values,
    labels=labels_plot,
    autopct="%1.1f%%",
    colors=wedge_colors,
    startangle=90,
    pctdistance=0.75,
    textprops={"fontsize": 9},
)
for at in autotexts:
    at.set_fontsize(9)
axes[1].set_title("Annotation source\nbreakdown (n = 268,034 cells)", fontsize=11)

fig.tight_layout()
save(fig, "Fig2C_annotation_validation.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 3 — Spatial domains (3A Novae, 3B Banksy, 3C Consensus)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 3 — Spatial domain maps...")

coords  = pd.read_csv(os.path.join(SPATIAL, "F1_spatial_domains", "coords.csv"))
novae   = pd.read_csv(os.path.join(SPATIAL, "F1_spatial_domains", "novae_domains.csv"))
banksy  = pd.read_csv(os.path.join(SPATIAL, "F1_spatial_domains", "banksy_domains.csv"))
cons    = pd.read_csv(os.path.join(SPATIAL, "F1_spatial_domains", "consensus_domains.csv"))

# merge all on cell_id
dom = coords.merge(novae[["cell_id","novae_domain"]], on="cell_id")
dom = dom.merge(banksy[["cell_id","banksy_domain"]], on="cell_id")
dom = dom.merge(cons[["cell_id","consensus_domain","agreement"]], on="cell_id")

# Biologically meaningful labels for Novae domains (top 10)
novae_labels = {
    "D998":  "Alveolar/Epithelial",
    "D1013": "Immune infiltrate",
    "D1000": "Stromal",
    "D1010": "Tumor core",
    "D1008": "Vascular",
    "D1004": "Mixed immune",
    "D1012": "Transitional",
    "D971":  "Lymphoid",
    "D1006": "Fibrotic",
    "D985":  "Rare niche",
}
dom["novae_label"] = dom["novae_domain"].map(novae_labels).fillna(dom["novae_domain"])

# Banksy: use domain IDs directly (14 domains)
dom["banksy_label"] = dom["banksy_domain"]

DOMAIN_PALETTE_10 = [
    "#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261",
    "#264653","#A8DADC","#6D6875","#B5838D","#E5989B"
]
DOMAIN_PALETTE_14 = [
    "#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261",
    "#264653","#A8DADC","#6D6875","#B5838D","#E5989B",
    "#90E0EF","#023E8A","#80B918","#C77DFF"
]
DISAGREE_COLOR = "#D3D3D3"

def plot_domain_map(ax, x, y, labels, palette, title, panel_letter):
    unique = [l for l in pd.Series(labels).value_counts().index]  # ordered by count
    color_map = {lbl: palette[i % len(palette)] for i, lbl in enumerate(unique)}
    for lbl in reversed(unique):  # plot smaller clusters on top
        mask = labels == lbl
        ax.scatter(x[mask], y[mask], c=color_map[lbl],
                   s=0.2, alpha=0.7, linewidths=0, rasterized=True, label=lbl)
    handles = [mpatches.Patch(color=color_map[l], label=l) for l in unique]
    ax.legend(handles=handles, fontsize=6, loc="upper right",
              framealpha=0.9, ncol=1, markerscale=3, handlelength=1)
    ax.set_title(title, fontsize=11)
    ax.set_aspect("equal")
    ax.axis("off")
    label_panel(ax, panel_letter)

# 3A — Novae
fig, ax = plt.subplots(figsize=(8, 7))
plot_domain_map(ax, dom["x"].values, dom["y"].values,
                dom["novae_label"].values, DOMAIN_PALETTE_10,
                "Spatial domains — foundation model (Novae)", "A")
fig.tight_layout()
save(fig, "Fig3A_novae_domains.png")

# 3B — Banksy
fig, ax = plt.subplots(figsize=(8, 7))
plot_domain_map(ax, dom["x"].values, dom["y"].values,
                dom["banksy_label"].values, DOMAIN_PALETTE_14,
                "Spatial domains — graph-based clustering (Banksy)", "B")
fig.tight_layout()
save(fig, "Fig3B_banksy_domains.png")

# 3C — Consensus (agreement cells colored, disagreement gray)
fig, ax = plt.subplots(figsize=(8, 7))
agree_mask   = dom["agreement"].values
disagree_mask= ~agree_mask

# plot disagreement first (background)
ax.scatter(dom["x"].values[disagree_mask], dom["y"].values[disagree_mask],
           c=DISAGREE_COLOR, s=0.1, alpha=0.3, linewidths=0, rasterized=True)

cons_unique = [l for l in dom.loc[agree_mask, "consensus_domain"].value_counts().index
               if l != "Disagreement"]
cons_colors = {lbl: DOMAIN_PALETTE_10[i % 10] for i, lbl in enumerate(cons_unique)}
for lbl in reversed(cons_unique):
    mask = agree_mask & (dom["consensus_domain"] == lbl).values
    ax.scatter(dom["x"].values[mask], dom["y"].values[mask],
               c=cons_colors[lbl], s=0.25, alpha=0.85,
               linewidths=0, rasterized=True)

handles = [mpatches.Patch(color=cons_colors[l], label=novae_labels.get(l, l))
           for l in cons_unique]
handles.append(mpatches.Patch(color=DISAGREE_COLOR, alpha=0.5,
                               label="Biologically plastic\n(methods disagree)"))
ax.legend(handles=handles, fontsize=6.5, loc="upper right",
          framealpha=0.9, markerscale=3, handlelength=1)
n_agree = agree_mask.sum()
ax.set_title(f"Consensus spatial domains\n"
             f"({n_agree:,} cells, {n_agree/len(dom)*100:.1f}% agreement)", fontsize=11)
ax.set_aspect("equal")
ax.axis("off")
label_panel(ax, "C")
fig.tight_layout()
save(fig, "Fig3C_consensus_domains.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 4B — Validated CCC axes (bar chart with p_adj Spacia)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 4B — Spacia validated axes...")

spacia = pd.read_csv(os.path.join(RESULTS, "02_biology", "phase_f_spacia",
                                  "F3_spacia_results.csv"))
validated = spacia[spacia["spatially_validated"] == True].copy()
validated["lr_axis"] = validated["ligand_complex"] + " → " + validated["receptor_complex"]
validated["neg_log_padj"] = -np.log10(validated["pathway_pval_adj"].clip(lower=1e-80))
validated["sender_clean"] = validated["source"].str.replace("_", " ")

# Group by LR axis
axis_order = ["ADAM17 → MUC1", "LTF → AGER", "CDH1 → EGFR", "S100B → AGER"]
axis_colors = {"ADAM17 → MUC1": "#C62828",
               "LTF → AGER":    "#1565C0",
               "CDH1 → EGFR":   "#2E7D32",
               "S100B → AGER":  "#6A1B9A"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Panel left: scatter of all validated interactions
ax = axes[0]
for axis_name in axis_order:
    lg, rc = axis_name.split(" → ")
    sub = validated[(validated["ligand_complex"] == lg) &
                    (validated["receptor_complex"] == rc)]
    ax.scatter(sub["scaled_weight"], sub["neg_log_padj"],
               c=axis_colors[axis_name], s=60, alpha=0.85,
               label=axis_name, zorder=3, edgecolors="white", linewidths=0.5)
ax.set_xlabel("LIANA+ scaled interaction score", fontsize=11)
ax.set_ylabel("−log₁₀(p_adj) Spacia", fontsize=11)
ax.set_title("Dual-validated interactions\n(LIANA+ × Spacia Bayesian MIL)", fontsize=11)
ax.legend(fontsize=9, framealpha=0.9)
ax.spines[["top","right"]].set_visible(False)
ax.axhline(-np.log10(0.05), color="gray", lw=1, ls="--", alpha=0.6)
ax.text(ax.get_xlim()[0], -np.log10(0.05)+0.5, "p_adj = 0.05",
        color="gray", fontsize=8)
label_panel(ax, "A")  # this becomes 4A complement

# Panel right: lollipop per LR axis showing max significance + n senders
ax2 = axes[1]
axis_summary = []
for axis_name in axis_order:
    lg, rc = axis_name.split(" → ")
    sub = validated[(validated["ligand_complex"] == lg) &
                    (validated["receptor_complex"] == rc)]
    axis_summary.append({
        "axis": axis_name,
        "n_senders": len(sub),
        "max_neg_log": sub["neg_log_padj"].max(),
        "color": axis_colors[axis_name],
    })
summary_df = pd.DataFrame(axis_summary).sort_values("max_neg_log")

y_pos = range(len(summary_df))
ax2.barh(list(y_pos), summary_df["max_neg_log"],
         color=summary_df["color"].tolist(), height=0.55, alpha=0.85)
ax2.set_yticks(list(y_pos))
ax2.set_yticklabels(summary_df["axis"].tolist(), fontsize=11)
for i, (_, row) in enumerate(summary_df.iterrows()):
    ax2.text(row["max_neg_log"] + 0.5, i,
             f"n={int(row['n_senders'])} senders",
             va="center", fontsize=9.5, fontweight="bold")
ax2.set_xlabel("Max −log₁₀(p_adj) Spacia", fontsize=11)
ax2.set_title("Dominant signaling axes\n(strongest spatially validated interaction)", fontsize=11)
ax2.spines[["top","right"]].set_visible(False)
label_panel(ax2, "B")

fig.tight_layout()
save(fig, "Fig4B_spacia_axes.png")


# ═══════════════════════════════════════════════════════════════════════════════
# FIG 5B — Cell-type proportions along pseudotime (Lineage 1)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig 5B — Pseudotime cell-type proportions...")

from scipy.ndimage import uniform_filter1d

pt_df = pd.read_csv(os.path.join(SPATIAL, "F5_pseudotime",
                                  "pseudotime_with_spatial.csv"), index_col=0)
lin1 = pt_df.dropna(subset=["slingPseudotime_1"]).copy()
lin1["pt_bin"] = pd.cut(lin1["slingPseudotime_1"], bins=20, labels=False)

# proportions per bin
prop = (lin1.groupby(["pt_bin", "cell_type_L2"])
            .size()
            .unstack(fill_value=0))
prop_norm = prop.div(prop.sum(axis=1), axis=0)

# top 8 by mean abundance
top_cts = prop_norm.mean().nlargest(8).index.tolist()
prop_norm = prop_norm[top_cts]

bin_edges = pd.cut(lin1["slingPseudotime_1"], bins=20, retbins=True)[1]
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# populations used to build the analysis (from F5_summary.json)
all_pops_pt = [
    "Tumor_proliferating", "Tumor_resting",
    "Macrophage_M1", "Macrophage_M2", "Monocyte_classical",
    "CD8_T_cytotoxic", "CD8_T_exhausted", "CD4_T_helper", "Treg", "NK_cytotoxic",
    "cDC1", "cDC2",
    "Endothelial_blood", "Endothelial_lymphatic",
    "Epithelial_general", "AT1", "Ciliated",
]

ct_colors_pt = {ct: ct_color.get(ct, "#999999") for ct in top_cts}
smoothed = {}
for ct in top_cts:
    y = prop_norm[ct].reindex(range(20), fill_value=0).values
    smoothed[ct] = uniform_filter1d(y.astype(float), size=3)
dominant_per_bin = pd.DataFrame(smoothed, index=bin_centers).idxmax(axis=1)

# right margin for end-of-line labels
fig, ax = plt.subplots(figsize=(13, 6))
fig.subplots_adjust(left=0.07, right=0.70, top=0.85, bottom=0.22)

for ct in top_cts:
    y_s = smoothed[ct]
    ax.plot(bin_centers, y_s, lw=2.5, color=ct_colors_pt[ct], zorder=3)
    mask = (dominant_per_bin == ct).values
    ax.fill_between(bin_centers, y_s, where=mask, alpha=0.18,
                    color=ct_colors_pt[ct], zorder=2)

# end-of-line labels spread vertically (no overlap)
xmax_val = bin_centers[-1]
end_vals = sorted([(smoothed[ct][-1], ct) for ct in top_cts])
label_positions = []
for ev, ct in end_vals:
    pos = ev
    for prev in label_positions:
        if abs(pos - prev) < 0.045:
            pos = prev + 0.045
    label_positions.append(pos)

for (ev, ct), lpos in zip(end_vals, label_positions):
    col = ct_colors_pt[ct]
    ax.annotate(ct.replace("_", " "),
                xy=(xmax_val, ev), xytext=(xmax_val * 1.04, lpos),
                xycoords="data", textcoords="data",
                arrowprops=dict(arrowstyle="-", color=col, lw=0.7, alpha=0.6),
                ha="left", va="center", fontsize=8.5, color=col, fontweight="bold",
                annotation_clip=False)

# trajectory arrow + start/end labels below x-axis
ax.annotate("", xy=(1.0, -0.22), xytext=(0.0, -0.22),
            xycoords="axes fraction", textcoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", color="#444444", lw=1.5))
ax.text(0.0, -0.26, "Tumor proliferating  (start)", transform=ax.transAxes,
        ha="left", va="top", fontsize=8.5, color="#e41a1c", fontweight="bold")
ax.text(1.0, -0.26, "Treg  (Lineage 1 end)", transform=ax.transAxes,
        ha="right", va="top", fontsize=8.5, color="#6a3d9a", fontweight="bold")

ax.axvline(bin_centers[0],  color="#e41a1c", lw=1.0, ls="--", alpha=0.4, zorder=1)
ax.axvline(bin_centers[-1], color="#6a3d9a", lw=1.0, ls="--", alpha=0.4, zorder=1)

# populations included box below trajectory arrow
rows_pt = [all_pops_pt[i:i+9] for i in range(0, len(all_pops_pt), 9)]
pop_lines = ["Included (17 cell types):  " + "   ".join(p.replace("_", " ")
             for p in rows_pt[0])]
for r in rows_pt[1:]:
    pop_lines.append("                              " + "   ".join(
        p.replace("_", " ") for p in r))
ax.text(0.0, -0.42, "\n".join(pop_lines), transform=ax.transAxes,
        ha="left", va="top", fontsize=7.5, color="#444444", linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#cccccc", lw=0.7))

ax.set_xlabel("Pseudotime (tumor → stroma gradient)", fontsize=11, labelpad=8)
ax.set_ylabel("Cell-type proportion", fontsize=11)
ax.set_title(
    f"Cell-type composition along spatial pseudotime — Lineage 1\n"
    f"(n = {len(lin1):,} cells in this lineage  |  {len(pt_df):,} total cells"
    f"  |  8 lineages total)",
    fontsize=11, pad=10)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xlim(bin_centers[0], bin_centers[-1])
ax.set_ylim(bottom=0)
label_panel(ax, "B")
fig.savefig(os.path.join(FIGDIR, "Fig5B_pseudotime_proportions.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# FIG S1 — QC metrics
# ═══════════════════════════════════════════════════════════════════════════════
print("Fig S1 — QC metrics...")

fig, axes = plt.subplots(1, 3, figsize=(13, 4))

def plot_hist(ax, data, xlabel, title, color, vline=None):
    ax.hist(data, bins=60, color=color, alpha=0.8, edgecolor="white", linewidth=0.3)
    if vline is not None:
        ax.axvline(vline, color="crimson", lw=1.5, ls="--",
                   label=f"median = {vline:.0f}")
        ax.legend(fontsize=9)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Number of cells", fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.spines[["top","right"]].set_visible(False)

tc = obs["transcript_counts"].values
ng = obs["n_genes_by_counts"].values
ca = obs["cell_area"].values

plot_hist(axes[0], tc, "Transcripts per cell", "Transcripts per cell",
          "#42A5F5", vline=np.median(tc))
plot_hist(axes[1], ng, "Genes detected per cell", "Genes detected per cell",
          "#66BB6A", vline=np.median(ng))
plot_hist(axes[2], ca, "Cell area (pixels²)", "Cell area",
          "#FFA726", vline=np.median(ca))

for i, ax in enumerate(axes):
    label_panel(ax, ["A","B","C"][i])

fig.suptitle(f"Quality metrics — {len(obs):,} cells after filtering", fontsize=13)
fig.tight_layout()
save(fig, "FigS1_qc_metrics.png")


print("\n✓ All manuscript figures generated.")
print(f"Output directory: {FIG_OUT}")
