#!/usr/bin/env python
"""M1 — Pseudobulk consensus across decoupler-py + muscat (edgeR/DESeq2/limma).

Reads the four DE tables produced upstream:

    pseudobulk_decoupler/<organ>/<celltype>/<contrast>__decoupler.tsv
    pseudobulk_muscat/<organ>/<celltype>/<contrast>__edgeR.tsv
    pseudobulk_muscat/<organ>/<celltype>/<contrast>__DESeq2.tsv
    pseudobulk_muscat/<organ>/<celltype>/<contrast>__limma_voom.tsv

For every (organ, celltype, contrast, gene), produces:

    consensus/<organ>/<celltype>/<contrast>__consensus.tsv
        gene, n_methods_called, methods_called, methods_failed,
        mean_log2FC, median_log2FC, log2FC_concordance,
        max_padj, min_padj, fisher_combined_p, fdr_combined,
        consensus_call ∈ {strict_consensus, majority, single, discordant}

    summary/method_concordance.tsv  — per-method pairwise concordance
    summary/venn_<organ>_<celltype>_<contrast>.png  — 4-set Venn

Consensus rules:
    strict_consensus : called by all 4 methods, log2FC sign matches in all 4
    majority         : called by ≥3 methods, log2FC sign matches
    single           : called by exactly 1 method
    discordant       : called by ≥2 methods but log2FC signs disagree

This script does NOT redo statistics. It only joins, harmonizes, and
reports concordance — purpose is the discussion section of the manuscript
("alcances, diferencias y consenso al que se puede llegar con cada
herramienta").

CLI
---
    python pseudobulk_consensus.py
        --decoupler-dir TBDs/cohort/results/pseudobulk_decoupler
        --muscat-dir    TBDs/cohort/results/pseudobulk_muscat
        --outdir        TBDs/cohort/results/pseudobulk_consensus
        [--padj-threshold 0.05] [--lfc-threshold 0.585]
        [--skip-venn]
"""
from __future__ import annotations

import argparse
import logging
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

log = logging.getLogger("pseudobulk_consensus")

METHODS = ("decoupler", "edgeR", "DESeq2", "limma_voom")


# ────────────────────────────────────────────────────────────────────────
def _load_method_tables(
    decoupler_dir: Path, muscat_dir: Path
) -> pd.DataFrame:
    """Stack all DE tables (long format) with a `method` column."""
    frames: List[pd.DataFrame] = []
    for path in decoupler_dir.rglob("*__decoupler.tsv"):
        df = pd.read_csv(path, sep="\t")
        df["method"] = "decoupler"
        organ = path.parents[1].name
        df["organ"] = organ
        df["celltype"] = path.parent.name
        frames.append(df)
    for method in ("edgeR", "DESeq2", "limma_voom"):
        for path in muscat_dir.rglob(f"*__{method}.tsv"):
            df = pd.read_csv(path, sep="\t")
            df["method"] = method
            organ = path.parents[1].name
            df["organ"] = organ
            df["celltype"] = path.parent.name
            frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"no DE tables under {decoupler_dir} or {muscat_dir}"
        )
    df = pd.concat(frames, ignore_index=True)
    needed = {"gene", "log2FoldChange", "padj", "contrast",
              "method", "organ", "celltype"}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(f"DE tables missing columns: {missing}")
    return df


# ────────────────────────────────────────────────────────────────────────
def _consensus_call(row: pd.Series, padj_threshold: float,
                    lfc_threshold: float) -> str:
    """Classify a single gene's cross-method behavior."""
    called: List[str] = []
    signs: List[int] = []
    for m in METHODS:
        p = row.get(f"padj_{m}")
        lfc = row.get(f"log2FC_{m}")
        if pd.notna(p) and pd.notna(lfc) and p < padj_threshold \
                and abs(lfc) >= lfc_threshold:
            called.append(m)
            signs.append(int(np.sign(lfc)))
    if not called:
        return "not_called"
    if len(called) == 4 and len(set(signs)) == 1:
        return "strict_consensus"
    if len(called) >= 3 and len(set(signs)) == 1:
        return "majority"
    if len(called) == 1:
        return "single"
    if len(set(signs)) > 1:
        return "discordant"
    return "partial"


def consensus_table(
    df_long: pd.DataFrame,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 0.585,
) -> pd.DataFrame:
    """Pivot per-method and compute consensus calls + Fisher-combined p."""
    # Pivot: one row per (organ, celltype, contrast, gene), columns per method
    pivot = df_long.pivot_table(
        index=["organ", "celltype", "contrast", "gene"],
        columns="method",
        values=["log2FoldChange", "padj"],
        aggfunc="first",
    )
    pivot.columns = [f"{val}_{m}" for val, m in pivot.columns]
    pivot = pivot.reset_index().rename(columns={
        f"log2FoldChange_{m}": f"log2FC_{m}" for m in METHODS
    })

    # Cross-method summaries
    lfc_cols = [c for c in pivot.columns if c.startswith("log2FC_")]
    padj_cols = [c for c in pivot.columns if c.startswith("padj_")]
    pivot["mean_log2FC"] = pivot[lfc_cols].mean(axis=1, skipna=True)
    pivot["median_log2FC"] = pivot[lfc_cols].median(axis=1, skipna=True)
    pivot["log2FC_sign_concordance"] = pivot[lfc_cols].apply(
        lambda r: 1.0 if len(set(np.sign(r.dropna()))) <= 1 else 0.0, axis=1
    )
    pivot["max_padj"] = pivot[padj_cols].max(axis=1, skipna=True)
    pivot["min_padj"] = pivot[padj_cols].min(axis=1, skipna=True)
    pivot["n_methods_tested"] = pivot[padj_cols].notna().sum(axis=1)

    # Fisher combined p-value across methods
    def fisher_combine(row: pd.Series) -> float:
        ps = [row[c] for c in padj_cols if pd.notna(row[c]) and row[c] > 0]
        if not ps:
            return np.nan
        return st.combine_pvalues(ps, method="fisher")[1]

    pivot["fisher_combined_p"] = pivot.apply(fisher_combine, axis=1)
    mask = pivot["fisher_combined_p"].notna()
    pivot.loc[mask, "fisher_combined_fdr"] = multipletests(
        pivot.loc[mask, "fisher_combined_p"], method="fdr_bh"
    )[1]

    pivot["consensus_call"] = pivot.apply(
        _consensus_call, axis=1,
        padj_threshold=padj_threshold, lfc_threshold=lfc_threshold,
    )
    return pivot


# ────────────────────────────────────────────────────────────────────────
def pairwise_method_concordance(
    df_long: pd.DataFrame,
    padj_threshold: float = 0.05,
) -> pd.DataFrame:
    """For each pair of methods, % of shared significant genes (Jaccard)."""
    rows: List[Dict] = []
    sig = df_long[df_long["padj"] < padj_threshold].copy()
    for (organ, ct, contrast), grp in sig.groupby(
            ["organ", "celltype", "contrast"]):
        gene_sets = {
            m: set(grp.loc[grp["method"] == m, "gene"]) for m in METHODS
        }
        for m1, m2 in combinations(METHODS, 2):
            a, b = gene_sets[m1], gene_sets[m2]
            jaccard = (len(a & b) / len(a | b)) if (a | b) else np.nan
            rows.append({
                "organ": organ, "celltype": ct, "contrast": contrast,
                "method_a": m1, "method_b": m2,
                "n_a": len(a), "n_b": len(b),
                "n_intersection": len(a & b),
                "n_union": len(a | b),
                "jaccard": jaccard,
            })
    return pd.DataFrame(rows)


# ────────────────────────────────────────────────────────────────────────
def venn4_png(
    sets: Dict[str, set], title: str, out: Path
) -> Optional[Path]:
    try:
        from matplotlib_venn import venn3                # type: ignore
        import matplotlib.pyplot as plt                  # type: ignore
    except ImportError:
        log.warning("matplotlib_venn not installed — skipping venn for %s", title)
        return None
    # 4-way venn isn't supported in matplotlib_venn — collapse limma_voom into edgeR for the figure
    fig, ax = plt.subplots(figsize=(5, 4))
    venn3(
        [sets["decoupler"], sets["edgeR"] | sets["limma_voom"], sets["DESeq2"]],
        ("decoupler", "edgeR+limma", "DESeq2"),
        ax=ax,
    )
    ax.set_title(title)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ────────────────────────────────────────────────────────────────────────
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--decoupler-dir", type=Path, required=True)
    p.add_argument("--muscat-dir", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--padj-threshold", type=float, default=0.05)
    p.add_argument("--lfc-threshold", type=float, default=0.585)
    p.add_argument("--skip-venn", action="store_true")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    df_long = _load_method_tables(args.decoupler_dir, args.muscat_dir)
    log.info("loaded %d rows across %d methods", len(df_long),
             df_long["method"].nunique())

    cons = consensus_table(
        df_long,
        padj_threshold=args.padj_threshold,
        lfc_threshold=args.lfc_threshold,
    )
    cons_dir = args.outdir / "consensus"
    cons_dir.mkdir(parents=True, exist_ok=True)
    for (organ, ct, contrast), grp in cons.groupby(
            ["organ", "celltype", "contrast"]):
        d = cons_dir / organ / ct
        d.mkdir(parents=True, exist_ok=True)
        out = d / f"{contrast}__consensus.tsv"
        grp.to_csv(out, sep="\t", index=False)
    log.info("wrote consensus tables under %s", cons_dir)

    summary_dir = args.outdir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    concord = pairwise_method_concordance(df_long, args.padj_threshold)
    concord.to_csv(summary_dir / "method_concordance.tsv",
                   sep="\t", index=False)
    log.info("wrote %s (%d rows)",
             summary_dir / "method_concordance.tsv", len(concord))

    if not args.skip_venn:
        for (organ, ct, contrast), grp in df_long.groupby(
                ["organ", "celltype", "contrast"]):
            sig = grp[grp["padj"] < args.padj_threshold]
            sets = {m: set(sig.loc[sig["method"] == m, "gene"]) for m in METHODS}
            if all(len(s) == 0 for s in sets.values()):
                continue
            out = summary_dir / f"venn_{organ}_{ct}_{contrast}.png"
            venn4_png(sets, f"{organ} | {ct} | {contrast}", out)

    log.info("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
