#!/usr/bin/env python3
"""
F4 — NicheDE Visualizations

Generates publication-quality figures from NicheDE results:
  1. Bubble plot: cell type × domain (n significant genes, color = mean padj)
  2. Top niche-DE genes heatmap per cell type
  3. Spatial map of top niche-DE gene expression per domain
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import spatialdata as sd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

# Crameri-inspired CVD-friendly palette for cell types
CT_COLORS = {
    "B_cell":      "#4C8BB5",
    "Endothelial": "#E07B54",
    "Epithelial":  "#5BA85B",
    "Macrophage":  "#9B6BBF",
    "NK_cell":     "#E8C53A",
    "T_cell":      "#E85A5A",
    "Tumor":       "#4ABFBF",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--results_dir", required=True, type=Path)
    p.add_argument("--sdata",       required=True, type=Path)
    p.add_argument("--outdir",      required=True, type=Path)
    return p.parse_args()


def plot_bubble(df: pd.DataFrame, outpath: Path) -> None:
    """Bubble plot: cell_type (x) × niche/domain (y), bubble = n_sig_genes."""
    pivot = df.groupby(["cell_type", "niche"]).agg(
        n_sig=("padj", lambda x: (x < 0.05).sum()),
        mean_padj=("padj", "mean")
    ).reset_index()

    pivot = pivot[pivot["n_sig"] > 0]
    if pivot.empty:
        log.warning("No significant niche-DE genes to plot.")
        return

    ct_order = sorted(pivot["cell_type"].unique())
    niche_order = sorted(pivot["niche"].unique())

    fig, ax = plt.subplots(figsize=(max(6, len(ct_order) * 1.2), max(5, len(niche_order) * 0.6)))

    for _, row in pivot.iterrows():
        x = ct_order.index(row["cell_type"])
        y = niche_order.index(row["niche"])
        size = np.sqrt(row["n_sig"]) * 60
        color = CT_COLORS.get(row["cell_type"], "#888888")
        alpha = 1 - row["mean_padj"]
        ax.scatter(x, y, s=size, c=color, alpha=max(0.3, min(1.0, alpha)),
                   edgecolors="white", linewidths=0.5, zorder=3)
        if row["n_sig"] >= 3:
            ax.text(x, y, str(int(row["n_sig"])), ha="center", va="center",
                    fontsize=6, color="white", fontweight="bold")

    ax.set_xticks(range(len(ct_order)))
    ax.set_xticklabels(ct_order, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(niche_order)))
    ax.set_yticklabels(niche_order, fontsize=9)
    ax.set_xlabel("Cell Type", fontsize=11)
    ax.set_ylabel("Spatial Domain (Banksy)", fontsize=11)
    ax.set_title("Niche-DE: Significant Genes per Cell Type × Domain\n"
                 "(bubble size = n genes, padj < 0.05)", fontsize=11)
    ax.grid(alpha=0.25, linestyle="--")
    ax.set_xlim(-0.7, len(ct_order) - 0.3)
    ax.set_ylim(-0.7, len(niche_order) - 0.3)

    # Legend for bubble size
    for n_ref in [1, 5, 20]:
        ax.scatter([], [], s=np.sqrt(n_ref)*60, c="gray", alpha=0.6,
                   label=f"n={n_ref}", edgecolors="white")
    ax.legend(title="Sig. genes", loc="upper right", fontsize=8, title_fontsize=8)

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Bubble plot: {outpath}")


def plot_top_genes_heatmap(df: pd.DataFrame, outpath: Path, top_n: int = 5) -> None:
    """Heatmap: top N niche-DE genes per (cell_type × domain) combination."""
    sig = df[df["padj"] < 0.05].copy()
    if sig.empty:
        log.warning("No padj<0.05 genes for heatmap.")
        return

    # Top genes per cell_type
    top = (sig.sort_values("padj")
             .groupby("cell_type")
             .head(top_n)
             .drop_duplicates("gene"))

    ct_order = sorted(top["cell_type"].unique())
    genes = top["gene"].unique()

    # Pivot: gene × cell_type (value = -log10 padj)
    top["log_padj"] = -np.log10(top["padj"].clip(1e-300))
    pivot = top.pivot_table(index="gene", columns="cell_type", values="log_padj",
                            aggfunc="max", fill_value=0)
    pivot = pivot.reindex(columns=[c for c in ct_order if c in pivot.columns])
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(max(6, len(ct_order) * 1.2), max(5, len(pivot) * 0.3)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xlabel("Cell Type", fontsize=11)
    ax.set_title(f"Top Niche-DE Genes per Cell Type\n(-log10 padj)", fontsize=11)
    plt.colorbar(im, ax=ax, shrink=0.6, label="-log10(padj)")
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Gene heatmap: {outpath}")


def plot_spatial_niche_map(sig_df: pd.DataFrame, sdata_path: Path, outpath: Path,
                           top_ct: str = "Macrophage") -> None:
    """Spatial map coloring cells by their niche-DE activity for a focal cell type."""
    ct_genes = sig_df[sig_df["cell_type"] == top_ct].sort_values("padj")
    if ct_genes.empty:
        log.warning(f"No niche-DE genes for {top_ct}; skipping spatial map.")
        return

    top_gene = ct_genes.iloc[0]["gene"]
    top_niche = ct_genes.iloc[0]["niche"]
    log.info(f"  Spatial map: {top_ct} / {top_gene} (niche={top_niche})")

    sdata = sd.read_zarr(str(sdata_path))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]

    coords = adata.obsm["spatial"]
    x, y = coords[:, 0], coords[:, 1]

    # Get expression of top gene
    gene_idx = list(adata.var_names).index(top_gene) if top_gene in adata.var_names else None
    if gene_idx is None:
        log.warning(f"Gene {top_gene} not found in adata.")
        return

    import scipy.sparse as sp_
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    expr = np.asarray(X[:, gene_idx].todense()).ravel() if sp_.issparse(X) else X[:, gene_idx]

    # Cell type mask
    ct_mask = adata.obs["cell_type"].values == top_ct

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: all cells colored by expression
    sc0 = axes[0].scatter(x, y, s=0.3, c=expr, cmap="Blues",
                           norm=mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=expr.max()),
                           rasterized=True)
    plt.colorbar(sc0, ax=axes[0], shrink=0.7, label="Raw counts")
    axes[0].set_title(f"{top_gene} expression (all cells)", fontsize=10)
    axes[0].set_aspect("equal")
    axes[0].axis("off")

    # Right: highlight focal cell type
    axes[1].scatter(x[~ct_mask], y[~ct_mask], s=0.2, c="lightgrey",
                    rasterized=True, label="Other")
    sc1 = axes[1].scatter(x[ct_mask], y[ct_mask], s=1.0,
                          c=expr[ct_mask], cmap="Reds",
                          norm=mcolors.PowerNorm(gamma=0.4, vmin=0,
                                                  vmax=max(expr[ct_mask].max(), 1)),
                          rasterized=True, label=top_ct)
    plt.colorbar(sc1, ax=axes[1], shrink=0.7, label="Raw counts")
    axes[1].set_title(f"{top_ct} — {top_gene}\n(top niche-DE gene, niche={top_niche})",
                      fontsize=10)
    axes[1].set_aspect("equal")
    axes[1].axis("off")

    fig.suptitle(f"Spatial map: top niche-DE gene for {top_ct}", fontsize=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Spatial map: {outpath}")


def _plot_no_results_summary(results_dir: Path, outdir: Path) -> None:
    """Informative summary figure when no significant niche-DE genes are found."""
    import json
    summary_path = results_dir / "nichede_summary.json"
    summary = {}
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)

    ct_tab_path = results_dir / "domain_celltype_crosstab.csv"
    has_tab = ct_tab_path.exists()

    fig, axes = plt.subplots(1, 2 if has_tab else 1,
                             figsize=(14 if has_tab else 7, 5))
    if not has_tab:
        axes = [axes]

    # Left panel: summary text
    ax = axes[0]
    ax.axis("off")
    lines = [
        "NicheDE — Phase D Results",
        "",
        f"Dataset: {summary.get('n_cells', '?'):,} cells × {summary.get('n_genes', '?')} genes",
        f"Cell types: {summary.get('n_cell_types', '?')}",
        f"σ (bandwidth): {summary.get('sigma', '?')} µm",
        f"Min cells per niche (C): {summary.get('C', '?')}",
        "",
        f"Genes with T-statistics: {summary.get('genes_with_tstats', 'N/A')}",
        f"Significant genes (p<0.05): {summary.get('significant_genes_p05', 0)}",
        "",
        "Interpretation:",
        "No significant niche-DE genes detected.",
        "This is expected for a 289-gene targeted panel.",
        "Niche-DE has maximum power on 5k+ gene panels",
        "(e.g., Xenium Prime) where subtle context-dependent",
        "expression differences are captured.",
    ]
    ax.text(0.05, 0.95, "\n".join(lines), transform=ax.transAxes,
            fontsize=10, va="top", fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="#F5F5F5", alpha=0.8))
    ax.set_title("NicheDE: No Significant Results", fontsize=12, fontweight="bold")

    # Right panel: domain × cell type heatmap
    if has_tab:
        ct_tab = pd.read_csv(ct_tab_path, index_col=0)
        ax2 = axes[1]
        im = ax2.imshow(ct_tab.values, aspect="auto", cmap="Blues",
                        interpolation="nearest")
        ax2.set_xticks(range(len(ct_tab.columns)))
        ax2.set_xticklabels(ct_tab.columns, rotation=35, ha="right", fontsize=8)
        ax2.set_yticks(range(len(ct_tab.index)))
        ax2.set_yticklabels(ct_tab.index, fontsize=8)
        ax2.set_xlabel("Cell Type", fontsize=10)
        ax2.set_ylabel("Banksy Domain", fontsize=10)
        ax2.set_title("Cell Type × Domain Distribution\n(subsample 50k cells)", fontsize=10)
        plt.colorbar(im, ax=ax2, shrink=0.7, label="Cell count")

    fig.suptitle("NicheDE Analysis — Xenium v1 289-Gene Panel", fontsize=13, fontweight="bold")
    fig.tight_layout()
    outpath = outdir / "nichede_summary_no_results.png"
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Summary figure: {outpath}")


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 70)
    log.info("F4 — NicheDE Visualizations")
    log.info("=" * 70)

    sig_path = args.results_dir / "nichede_significant_genes.csv"
    if not sig_path.exists():
        log.error(f"Results not found: {sig_path}")
        sys.exit(1)

    sig_df = pd.read_csv(sig_path)
    if sig_df.empty:
        log.warning("No significant niche-DE genes — generating summary figure.")
        _plot_no_results_summary(args.results_dir, args.outdir)
        return

    log.info(f"  Total niche-DE entries: {len(sig_df)}")
    log.info(f"  padj<0.05 genes: {(sig_df['padj'] < 0.05).sum()}")

    # 1. Bubble plot
    plot_bubble(sig_df, args.outdir / "nichede_bubble_plot.png")

    # 2. Heatmap
    plot_top_genes_heatmap(sig_df, args.outdir / "nichede_gene_heatmap.png", top_n=5)

    # 3. Spatial map — focal cell type with most niche-DE genes
    if "cell_type" in sig_df.columns and "padj" in sig_df.columns:
        ct_counts = sig_df[sig_df["padj"] < 0.05].groupby("cell_type").size()
        if not ct_counts.empty:
            focal_ct = ct_counts.idxmax()
            log.info(f"  Focal cell type for spatial map: {focal_ct} ({ct_counts[focal_ct]} genes)")
            plot_spatial_niche_map(
                sig_df, args.sdata,
                args.outdir / f"nichede_spatial_map_{focal_ct}.png",
                top_ct=focal_ct
            )

    log.info(f"\n✓ Visualizations complete → {args.outdir}")


if __name__ == "__main__":
    main()
