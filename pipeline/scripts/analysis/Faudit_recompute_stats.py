#!/usr/bin/env python3
"""
Faudit_recompute_stats.py — Audit fix script

Recomputes:
  1. LIANA+ BH FDR (correct multiple testing for CCC)
  2. Spacia Bonferroni (correct testing across 25 jobs)

Run from xenium_pipeline env:
  conda run -n xenium_pipeline python Faudit_recompute_stats.py
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.paths import project_root

BASE = project_root()
LIANA_ALL = BASE / "human_lung_cancer/results/02_biology/ccc_liana_L2_v3/liana_L2_all_interactions.csv"
SPACIA_RES = BASE / "human_lung_cancer/results/02_biology/phase_f_spacia/F3_spacia_results.csv"
OUT = BASE / "human_lung_cancer/results/02_biology/audit_corrected"
OUT.mkdir(parents=True, exist_ok=True)


def fix_liana_bh():
    print("=" * 60)
    print("LIANA+ — BH FDR correction")
    print("=" * 60)
    df = pd.read_csv(LIANA_ALL)
    n_total = len(df)

    # Exclude self-interactions (consistent with original script)
    df = df[df["source"] != df["target"]].copy()
    n_no_self = len(df)
    print(f"  Total interactions:        {n_total:,}")
    print(f"  After self-removal:        {n_no_self:,}")

    # Apply BH on cellphone_pvals
    pvals = df["cellphone_pvals"].fillna(1.0).values
    rej, padj, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    df["cellphone_pvals_BH_FDR"] = padj

    sig_strict = (padj < 0.05).sum()
    sig_lenient = (padj < 0.10).sum()

    # Compare with original "significant" count (top-rank-based, not FDR)
    print(f"\n  Significant @ BH FDR < 0.05:  {sig_strict:,}")
    print(f"  Significant @ BH FDR < 0.10:  {sig_lenient:,}")
    print(f"  Original (rank-based):        1,565  (NOT actually FDR)")

    df_sig = df[padj < 0.05].copy()
    out_csv = OUT / "liana_L2_significant_BH_FDR_corrected.csv"
    df_sig.to_csv(out_csv, index=False)
    print(f"\n  Saved: {out_csv}")

    # Per-celltype counts
    pair_counts = df_sig.groupby(["source", "target"]).size().rename("n").reset_index()
    pair_counts.to_csv(OUT / "liana_L2_pair_counts_BH_FDR.csv", index=False)

    return {
        "n_total": int(n_total),
        "n_no_self": int(n_no_self),
        "n_sig_BH_FDR_005": int(sig_strict),
        "n_sig_BH_FDR_010": int(sig_lenient),
    }


def fix_spacia_bonferroni():
    print("\n" + "=" * 60)
    print("Spacia — Bonferroni correction across all submitted tests")
    print("=" * 60)
    df = pd.read_csv(SPACIA_RES)
    n_submitted = len(df)
    df_run = df[df["pathway_pval"].notna()].copy()
    n_run = len(df_run)
    print(f"  Submitted:                 {n_submitted}")
    print(f"  Completed (had results):   {n_run}")
    print(f"  Failed (rare cell types):  {n_submitted - n_run}")

    # Original criterion: pathway_pval_adj < 0.05 (per-test, no cross-test correction)
    n_orig = (df_run["pathway_pval_adj"] < 0.05).sum()
    print(f"\n  Original (per-test pval_adj < 0.05): {n_orig}/{n_run}")

    # Bonferroni across N_RUN tests
    bonf_threshold_run = 0.05 / n_run
    bonf_threshold_submitted = 0.05 / n_submitted
    n_bonf_run = (df_run["pathway_pval"] < bonf_threshold_run).sum()
    n_bonf_submitted = (df_run["pathway_pval"] < bonf_threshold_submitted).sum()
    print(f"\n  Bonferroni across {n_run} completed (p<{bonf_threshold_run:.2e}):  {n_bonf_run}/{n_run}")
    print(f"  Bonferroni across {n_submitted} submitted (p<{bonf_threshold_submitted:.2e}): {n_bonf_submitted}/{n_run} of completed = {n_bonf_submitted}/{n_submitted} of attempted")

    # BH FDR for comparison
    pvals_run = df_run["pathway_pval"].values
    rej_bh, padj_bh, _, _ = multipletests(pvals_run, alpha=0.05, method="fdr_bh")
    df_run["pathway_pval_BH"] = padj_bh
    n_bh = rej_bh.sum()
    print(f"  BH FDR across {n_run} completed:    {n_bh}/{n_run}")

    # Save corrected table
    df_run["bonferroni_threshold_completed"] = bonf_threshold_run
    df_run["bonferroni_threshold_submitted"] = bonf_threshold_submitted
    df_run["validated_bonferroni_completed"] = df_run["pathway_pval"] < bonf_threshold_run
    df_run["validated_bonferroni_submitted"] = df_run["pathway_pval"] < bonf_threshold_submitted
    df_run["validated_BH_FDR_005"] = rej_bh

    out_csv = OUT / "spacia_results_corrected.csv"
    df_run.to_csv(out_csv, index=False)
    print(f"\n  Saved: {out_csv}")

    # Show what survives Bonferroni (most stringent)
    surviving = df_run[df_run["pathway_pval"] < bonf_threshold_submitted].sort_values("pathway_pval")
    print(f"\n  Survivors of Bonferroni (across {n_submitted} attempted, p<{bonf_threshold_submitted:.2e}):")
    print(surviving[["source", "target", "ligand_complex", "receptor_complex", "pathway_pval"]].to_string(index=False))

    return {
        "n_submitted": int(n_submitted),
        "n_completed": int(n_run),
        "n_failed": int(n_submitted - n_run),
        "n_orig_per_test_005": int(n_orig),
        "n_bonferroni_completed": int(n_bonf_run),
        "n_bonferroni_submitted": int(n_bonf_submitted),
        "n_BH_FDR_005": int(n_bh),
        "bonf_threshold_completed": float(bonf_threshold_run),
        "bonf_threshold_submitted": float(bonf_threshold_submitted),
    }


def main():
    summary = {
        "liana": fix_liana_bh(),
        "spacia": fix_spacia_bonferroni(),
    }
    with open(OUT / "audit_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n" + "=" * 60)
    print("AUDIT RECOMPUTE COMPLETE")
    print("=" * 60)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
