#!/usr/bin/env python3
"""
F2 — SVG Consensus: Multi-method Spatially Variable Gene integration

Integrates results from three independent SVG detection methods:
  1. Hotspot (autocorrelation, FDR-corrected)
  2. nnSVG   (LR statistic + padj)
  3. Moran's I (Phase 1A, already computed)

Consensus criterion: significant in ≥2 of 3 methods.
Reference: Salas et al. Nat Methods 22, 813-823 (2025).
"""

from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--hotspot_csv",  required=True, type=Path,
                   help="hotspot_svg_scores.csv from F2_hotspot_svg.py")
    p.add_argument("--nnsvg_csv",    required=True, type=Path,
                   help="nnsvg_svg_scores.csv from F2_nnsvg.R")
    p.add_argument("--morans_csv",   required=True, type=Path,
                   help="morans_i_all_genes.csv from Phase 1A")
    p.add_argument("--fdr_morans",   type=float, default=0.05,
                   help="FDR threshold for Moran's I significance (default: 0.05)")
    p.add_argument("--outdir",       required=True, type=Path)
    return p.parse_args()


def load_hotspot(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # columns: gene, C, Z, Pval, FDR, significant
    df = df[["gene", "C", "FDR", "significant"]].copy()
    df.columns = ["gene", "hotspot_C", "hotspot_FDR", "hotspot_sig"]
    return df.set_index("gene")


def load_nnsvg(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # columns: gene, rank, LR_stat, prop_sv, pval, padj, significant
    df = df[["gene", "rank", "LR_stat", "prop_sv", "padj", "significant"]].copy()
    df.columns = ["gene", "nnsvg_rank", "nnsvg_LR", "nnsvg_prop_sv", "nnsvg_padj", "nnsvg_sig"]
    return df.set_index("gene")


def load_morans(path: Path, fdr_threshold: float) -> pd.DataFrame:
    df = pd.read_csv(path)
    # columns: I, pval_norm, var_norm, pval_norm_fdr_bh, gene
    df = df[["gene", "I", "pval_norm_fdr_bh"]].copy()
    df.columns = ["gene", "morans_I", "morans_FDR"]
    df["morans_sig"] = df["morans_FDR"] < fdr_threshold
    return df.set_index("gene")


def plot_venn(sets: dict[str, set], outpath: Path) -> None:
    h_set = sets["Hotspot"]
    n_set = sets["nnSVG"]
    m_set = sets["Moran"]

    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn3([h_set, n_set, m_set], set_labels=("Hotspot", "nnSVG", "Moran's I"), ax=ax)
    ax.set_title("SVG Consensus — Method Overlap", fontsize=13, pad=12)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Venn saved: {outpath}")


def plot_consensus_heatmap(df: pd.DataFrame, outpath: Path, top_n: int = 50) -> None:
    """Dot/heatmap showing top consensus SVGs with scores from each method."""
    top = df.head(top_n).copy()
    # Normalize scores 0-1 per method for visualization
    for col, asc in [("hotspot_C", False), ("nnsvg_prop_sv", False), ("morans_I", False)]:
        if col in top.columns:
            mn, mx = top[col].min(), top[col].max()
            top[f"{col}_norm"] = (top[col] - mn) / (mx - mn + 1e-9)

    fig, axes = plt.subplots(1, 3, figsize=(13, max(6, top_n * 0.22)), sharey=True)
    methods = [
        ("hotspot_C_norm",    "Hotspot C",        "#E07B54"),
        ("nnsvg_prop_sv_norm","nnSVG prop_sv",     "#4C8BB5"),
        ("morans_I_norm",     "Moran's I",         "#5BA85B"),
    ]
    genes = top["gene"].tolist()
    y = np.arange(len(genes))

    for ax, (col, label, color) in zip(axes, methods):
        vals = top[col].values if col in top.columns else np.zeros(len(top))
        ax.barh(y, vals, color=color, alpha=0.75, height=0.7)
        ax.set_xlim(0, 1.05)
        ax.set_title(label, fontsize=10)
        ax.axvline(0, color="black", lw=0.5)
        if ax == axes[0]:
            ax.set_yticks(y)
            ax.set_yticklabels(genes, fontsize=7)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.3)

    fig.suptitle(f"Top {top_n} Consensus SVGs — Normalized Scores", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info(f"  Heatmap saved: {outpath}")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    log.info("=" * 70)
    log.info("F2 — SVG Multi-method Consensus")
    log.info("=" * 70)

    # ── Load results ──────────────────────────────────────────────────────────
    log.info("\n[1/4] Loading per-method results...")
    hs  = load_hotspot(args.hotspot_csv)
    nn  = load_nnsvg(args.nnsvg_csv)
    mi  = load_morans(args.morans_csv, args.fdr_morans)

    log.info(f"  Hotspot:  {len(hs):3d} genes | {hs['hotspot_sig'].sum():3d} significant")
    log.info(f"  nnSVG:    {len(nn):3d} genes | {nn['nnsvg_sig'].sum():3d} significant")
    log.info(f"  Moran's I:{len(mi):3d} genes | {mi['morans_sig'].sum():3d} significant")

    # ── Merge ─────────────────────────────────────────────────────────────────
    log.info("\n[2/4] Merging tables...")
    all_genes = hs.index.union(nn.index).union(mi.index)
    merged = pd.DataFrame(index=all_genes)
    merged = merged.join(hs, how="left").join(nn, how="left").join(mi, how="left")

    # Fill missing significances as False
    for col in ["hotspot_sig", "nnsvg_sig", "morans_sig"]:
        merged[col] = merged[col].fillna(False).astype(bool)

    # Count how many methods flag each gene
    merged["n_methods_significant"] = (
        merged["hotspot_sig"].astype(int) +
        merged["nnsvg_sig"].astype(int) +
        merged["morans_sig"].astype(int)
    )

    # Consensus: ≥2 of 3 methods
    merged["consensus_svg"] = merged["n_methods_significant"] >= 2
    merged["gene"] = merged.index

    n_consensus = int(merged["consensus_svg"].sum())
    log.info(f"  Consensus SVGs (≥2/3 methods): {n_consensus}/{len(merged)}")

    # ── Rank consensus genes ──────────────────────────────────────────────────
    # Sort by: n_methods DESC, then Hotspot C DESC (most robust first)
    consensus = merged[merged["consensus_svg"]].copy()
    consensus = consensus.sort_values(
        ["n_methods_significant", "hotspot_C"],
        ascending=[False, False]
    ).reset_index(drop=True)

    # ── Save outputs ──────────────────────────────────────────────────────────
    log.info("\n[3/4] Saving outputs...")

    # Main consensus table
    consensus_cols = [
        "gene", "n_methods_significant",
        "hotspot_C", "hotspot_FDR", "hotspot_sig",
        "nnsvg_rank", "nnsvg_LR", "nnsvg_prop_sv", "nnsvg_padj", "nnsvg_sig",
        "morans_I", "morans_FDR", "morans_sig",
    ]
    consensus_out = consensus[[c for c in consensus_cols if c in consensus.columns]]
    consensus_out.to_csv(args.outdir / "svg_consensus_table.csv", index=False)
    log.info(f"  Saved: svg_consensus_table.csv ({len(consensus_out)} genes)")

    # Full merged table (all genes, for supplementary)
    merged_cols = [
        "gene", "n_methods_significant", "consensus_svg",
        "hotspot_C", "hotspot_FDR", "hotspot_sig",
        "nnsvg_rank", "nnsvg_LR", "nnsvg_prop_sv", "nnsvg_padj", "nnsvg_sig",
        "morans_I", "morans_FDR", "morans_sig",
    ]
    merged_out = merged[[c for c in merged_cols if c in merged.columns]].reset_index(drop=True)
    merged_out.to_csv(args.outdir / "svg_all_genes.csv", index=False)
    log.info(f"  Saved: svg_all_genes.csv ({len(merged_out)} genes)")

    # Venn data JSON
    sets = {
        "Hotspot": set(hs.index[hs["hotspot_sig"]]),
        "nnSVG":   set(nn.index[nn["nnsvg_sig"]]),
        "Moran":   set(mi.index[mi["morans_sig"]]),
    }
    venn_data = {k: list(v) for k, v in sets.items()}
    venn_data["consensus"] = list(consensus["gene"])
    (args.outdir / "svg_venn_data.json").write_text(json.dumps(venn_data, indent=2))
    log.info("  Saved: svg_venn_data.json")

    # ── Figures ───────────────────────────────────────────────────────────────
    log.info("\n[4/4] Generating figures...")

    # Try Venn (requires matplotlib-venn)
    try:
        plot_venn(sets, args.outdir / "svg_venn.png")
    except ImportError:
        log.warning("  matplotlib-venn not installed; skipping Venn diagram")

    # Heatmap of top consensus SVGs
    plot_consensus_heatmap(consensus, args.outdir / "svg_consensus_heatmap.png", top_n=50)

    # ── Summary ───────────────────────────────────────────────────────────────
    log.info("\n  Top 20 Consensus SVGs (≥2/3 methods):")
    for _, row in consensus.head(20).iterrows():
        methods_str = "/".join([
            "HS" if row.get("hotspot_sig", False) else "--",
            "NN" if row.get("nnsvg_sig",   False) else "--",
            "MI" if row.get("morans_sig",   False) else "--",
        ])
        log.info(f"    {row['gene']:15s}  [{methods_str}]  "
                 f"C={row.get('hotspot_C', float('nan')):.3f}  "
                 f"I={row.get('morans_I',  float('nan')):.3f}")

    summary = {
        "n_genes_total":       int(len(merged)),
        "n_hotspot_sig":       int(hs["hotspot_sig"].sum()),
        "n_nnsvg_sig":         int(nn["nnsvg_sig"].sum()),
        "n_morans_sig":        int(mi["morans_sig"].sum()),
        "n_consensus_svg":     n_consensus,
        "n_all3_methods":      int((merged["n_methods_significant"] == 3).sum()),
        "n_exactly_2methods":  int((merged["n_methods_significant"] == 2).sum()),
        "execution_seconds":   round(time.time() - t0, 1),
    }
    (args.outdir / "svg_consensus_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    log.info(f"\n  Genes in all 3 methods : {summary['n_all3_methods']}")
    log.info(f"  Genes in exactly 2     : {summary['n_exactly_2methods']}")
    log.info(f"\n✓ SVG Consensus complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
