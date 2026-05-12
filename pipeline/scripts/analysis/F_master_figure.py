"""
Master Figure — Xenium Lung Cancer Pipeline
5-panel publication figure integrating all analytical phases.

Panel A: Spatial domains map (Novae + consensus agreement)
Panel B: SVG consensus heatmap (Hotspot/nnSVG/Moran's I)
Panel C: CCC bubble plot LIANA+ top interactions
Panel D: Pseudotime UMAP (Slingshot lineages)
Panel E: Annotation cross-validation ARI heatmap

Uses xenium_pipeline env.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from pathlib import Path
import logging
import sys

# ── Paths ─────────────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
RESULTS = BASE / "human_lung_cancer/results"
SPATIAL_DIR = RESULTS / "03_phase3_spatial"
FIG_DIR = RESULTS / "figures"
OUT_DIR = FIG_DIR  # flat structure post-2026-05-08 reorganization
OUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ── Colors ────────────────────────────────────────────────────────────────────
# Crameri-inspired palette for spatial domains (10 Novae domains)
DOMAIN_COLORS = [
    "#3B4CC0", "#6282EA", "#93B5E5", "#C7D8F0",
    "#F7F7F7", "#F5C7A5", "#E88061", "#C63F3F",
    "#7B2020", "#3D1C1C"
]

# Cell type colors (L2, 26 types) — publication-quality
CELLTYPE_COLORS = {
    "AT1": "#E63946", "Ciliated": "#F4A261", "Epithelial_general": "#E9C46A",
    "Smooth_muscle": "#2A9D8F", "Pericyte": "#264653",
    "Endothelial_blood": "#457B9D", "Endothelial_lymphatic": "#1D3557",
    "Tumor_resting": "#6D6875", "Tumor_proliferating": "#C77DFF",
    "CD4_T_helper": "#52B788", "CD8_T_cytotoxic": "#1B9E77",
    "CD8_T_exhausted": "#74C69D", "Treg": "#40916C",
    "NK_cytotoxic": "#B7E4C7", "B_naive": "#FFC8DD",
    "B_memory": "#FFAFCC", "Plasma": "#BDE0FE",
    "Plasmablast": "#A2D2FF", "Macrophage_M1": "#FB8500",
    "Macrophage_M2": "#FFB703", "Monocyte_classical": "#FD9E02",
    "cDC1": "#AE2012", "cDC2": "#BB3E03", "DC_mature": "#CA6702",
    "pDC": "#EE9B00", "Unassigned": "#CCCCCC",
}


def panel_A_spatial_domains(ax):
    """Spatial map colored by Novae domain, with consensus agreement highlighted."""
    log.info("Panel A: Spatial domains …")

    coords = pd.read_csv(SPATIAL_DIR / "F1_spatial_domains/coords.csv")
    domains = pd.read_csv(SPATIAL_DIR / "F1_spatial_domains/consensus_domains.csv")
    df = coords.merge(domains, on="cell_id")

    # Assign numeric domain index for coloring
    domain_labels = df["novae_domain"].unique()
    domain_map = {d: i for i, d in enumerate(sorted(domain_labels, key=str))}
    df["domain_idx"] = df["novae_domain"].map(domain_map)

    n_domains = len(domain_labels)
    sorted_labels = sorted(domain_labels, key=str)
    palette = (DOMAIN_COLORS * 3)[:n_domains]  # repeat if more domains than colors
    cmap = ListedColormap(palette)

    # Plot all cells, then overlay consensus cells
    pt_size = 0.3
    agreed = df[df["agreement"]]
    disagreed = df[~df["agreement"]]

    # Disagreement cells at low alpha (show domain structure, de-emphasized)
    ax.scatter(disagreed["x"], disagreed["y"],
               c=disagreed["domain_idx"], cmap=cmap, vmin=0, vmax=n_domains - 1,
               s=pt_size * 0.4, rasterized=True, linewidths=0, alpha=0.2)
    # Consensus cells at full opacity
    ax.scatter(agreed["x"], agreed["y"],
               c=agreed["domain_idx"], cmap=cmap, vmin=0, vmax=n_domains - 1,
               s=pt_size, rasterized=True, linewidths=0, alpha=0.9)

    legend_elements = [
        Patch(facecolor=palette[i], label=f"D{str(sorted_labels[i])[-3:]}")
        for i in range(n_domains)
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=5.5,
              framealpha=0.85, ncol=2, markerscale=2,
              title="Novae domains", title_fontsize=5.5)

    ax.set_aspect("equal")
    ax.axis("off")
    n_agreed = df["agreement"].sum()
    pct_agreed = 100 * n_agreed / len(df)
    ax.set_title(f"A  Spatial Domains (Novae 10 × Banksy 14)\n"
                 f"{n_agreed:,}/{len(df):,} strict-consensus ({pct_agreed:.0f}%) · "
                 f"ARI=0.15 NMI=0.31",
                 fontsize=9, fontweight="bold", loc="left")


def panel_B_svg_heatmap(ax):
    """SVG consensus heatmap from Phase C."""
    log.info("Panel B: SVG heatmap …")
    img_path = FIG_DIR / "03_svg_consensus_heatmap.png"
    if img_path.exists():
        img = mpimg.imread(str(img_path))
        ax.imshow(img, aspect="auto")
        ax.axis("off")
        ax.set_title("B  SVG Consensus\n(Hotspot · nnSVG · Moran's I)",
                     fontsize=9, fontweight="bold", loc="left")
    else:
        ax.text(0.5, 0.5, "SVG heatmap\n(not found)", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="gray")
        ax.axis("off")


def panel_C_ccc_bubble(ax):
    """CCC validation: Spacia cross-validation figure (if available), else LIANA+ bubble."""
    log.info("Panel C: CCC validation …")
    # Use the audited Bonferroni-corrected Spacia axes figure (Fig4B in manuscript)
    spacia_path = FIG_DIR / "ms_main_fig4B_spacia_axes.png"
    liana_path = FIG_DIR / "06_ccc_liana_L2_heatmap.png"
    img_path = spacia_path if spacia_path.exists() else liana_path
    if img_path.exists():
        img = mpimg.imread(str(img_path))
        ax.imshow(img, aspect="auto")
        ax.axis("off")
        ax.set_title("C  CCC Validation (LIANA+ × Spacia)\n"
                     "20/30 submitted Bonferroni-validated · ADAM17/MUC1 hub (6 senders)",
                     fontsize=9, fontweight="bold", loc="left")
    else:
        ax.text(0.5, 0.5, "CCC validation\n(not found)", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="gray")
        ax.axis("off")


def panel_D_pseudotime(ax):
    """Pseudotime panel — audited M1↔M2 macrophage continuum (post-2026-05-06)."""
    log.info("Panel D: Macrophage M1↔M2 pseudotime (audited) …")
    candidate_paths = [
        FIG_DIR / "07_pseudotime_paga_annotated.png",
        FIG_DIR / "07_pseudotime_concordance.png",
        FIG_DIR / "07_pseudotime_umap_overview.png",
    ]
    img_path = next((p for p in candidate_paths if p.exists()), None)
    if img_path is not None:
        img = mpimg.imread(str(img_path))
        ax.imshow(img, aspect="auto")
        ax.axis("off")
        ax.set_title("D  Macrophage M1↔M2 Pseudotime (Slingshot + DPT)\n"
                     "22,902 cells · Spearman r = 0.48 · 227 BH-FDR genes",
                     fontsize=9, fontweight="bold", loc="left")
    else:
        ax.text(0.5, 0.5, "Pseudotime UMAP\n(not found)", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="gray")
        ax.axis("off")


def panel_E_annotation_validation(ax):
    """Annotation cross-validation heatmap (ARI=0.742)."""
    log.info("Panel E: Annotation x-val …")
    img_path = FIG_DIR / "02_annot_ari_immune_granular.png"
    if img_path.exists():
        img = mpimg.imread(str(img_path))
        ax.imshow(img, aspect="auto")
        ax.axis("off")
        ax.set_title("E  Annotation Cross-Validation\n"
                     "L2 v3 vs immune subclustering · ARI = 0.742",
                     fontsize=9, fontweight="bold", loc="left")
    else:
        ax.text(0.5, 0.5, "Annotation x-val\n(not found)", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="gray")
        ax.axis("off")


def make_master_figure():
    """Assemble 5-panel master figure in a 2×3 grid (A large, B-E smaller)."""
    log.info("=" * 60)
    log.info("Building Master Figure")
    log.info("=" * 60)

    # Layout: A (left half, full height) | B top-right, C mid-right
    #         D bottom-left-half, E bottom-right-half
    fig = plt.figure(figsize=(20, 14), facecolor="white")

    # GridSpec: 2 rows × 3 cols
    # Row 0: A (spans cols 0-1), B (col 2 top), C (col 2 bottom top)
    # Row 1: D (cols 0-1), E (col 2)
    gs = gridspec.GridSpec(2, 3, figure=fig,
                           left=0.03, right=0.97,
                           top=0.93, bottom=0.03,
                           wspace=0.12, hspace=0.20)

    ax_A = fig.add_subplot(gs[0, 0:2])  # top-left 2/3
    ax_B = fig.add_subplot(gs[0, 2])    # top-right 1/3
    ax_C = fig.add_subplot(gs[1, 0])    # bottom-left 1/3
    ax_D = fig.add_subplot(gs[1, 1])    # bottom-mid 1/3
    ax_E = fig.add_subplot(gs[1, 2])    # bottom-right 1/3

    panel_A_spatial_domains(ax_A)
    panel_B_svg_heatmap(ax_B)
    panel_C_ccc_bubble(ax_C)
    panel_D_pseudotime(ax_D)
    panel_E_annotation_validation(ax_E)

    # Main title
    fig.suptitle(
        "Xenium Human Lung Cancer — Multi-Tool Spatial Transcriptomics Pipeline\n"
        "268,034 cells · 289 genes · 26 cell types (v3 hybrid annotation)",
        fontsize=12, fontweight="bold", y=0.98,
    )

    out_path = OUT_DIR / "master_figure.png"
    plt.savefig(str(out_path), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    log.info(f"Saved → {out_path}")

    # Also save PDF for submission
    out_pdf = OUT_DIR / "master_figure.pdf"
    fig2 = plt.figure(figsize=(20, 14), facecolor="white")
    gs2 = gridspec.GridSpec(2, 3, figure=fig2,
                            left=0.03, right=0.97,
                            top=0.93, bottom=0.03,
                            wspace=0.12, hspace=0.20)
    ax_A2 = fig2.add_subplot(gs2[0, 0:2])
    ax_B2 = fig2.add_subplot(gs2[0, 2])
    ax_C2 = fig2.add_subplot(gs2[1, 0])
    ax_D2 = fig2.add_subplot(gs2[1, 1])
    ax_E2 = fig2.add_subplot(gs2[1, 2])

    panel_A_spatial_domains(ax_A2)
    panel_B_svg_heatmap(ax_B2)
    panel_C_ccc_bubble(ax_C2)
    panel_D_pseudotime(ax_D2)
    panel_E_annotation_validation(ax_E2)

    fig2.suptitle(
        "Xenium Human Lung Cancer — Multi-Tool Spatial Transcriptomics Pipeline\n"
        "268,034 cells · 289 genes · 26 cell types (v3 hybrid annotation)",
        fontsize=12, fontweight="bold", y=0.98,
    )
    plt.savefig(str(out_pdf), bbox_inches="tight", facecolor="white")
    plt.close()
    log.info(f"Saved PDF → {out_pdf}")

    return str(out_path)


if __name__ == "__main__":
    make_master_figure()
