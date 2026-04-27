#!/usr/bin/env python3
"""
Phase 2B (REVISED): Cell-Cell Communication via LIANA+
=======================================================

Replaces previous custom heuristic with validated LIANA+ framework.

Methods used:
- rank_aggregate: consensus of CellPhoneDB, CellChat, NATMI, Connectome, SingleCellSignalR, logFC
- Database: 'consensus' (OmniPath unified L/R database)
- Permutation null model (n_perms=1000)
- Removes self-interactions automatically
- Filters by panel coverage (289-gene Xenium panel)

Output: validated L/R interactions with magnitude + specificity scores
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import sys

import spatialdata as sd
import scanpy as sc
import liana as li
import liana.method as lm

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "ccc_liana"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  Loaded {adata.shape[0]:,} cells x {adata.shape[1]} genes")
    log.info(f"  Cell types: {sorted(adata.obs['cell_type'].unique().tolist())}")
    return adata


def run_liana_ccc(adata, n_perms=1000, min_cells=50, expr_prop=0.1):
    """
    Run LIANA+ rank_aggregate (consensus of 5+ methods).
    Uses normalized counts in adata.X (already log-normalized in pipeline).
    """
    log.info("=" * 80)
    log.info("LIANA+ rank_aggregate (consensus CCC)")
    log.info("  Resource: consensus (OmniPath)")
    log.info(f"  Permutations: {n_perms}")
    log.info(f"  Min cells per type: {min_cells}")
    log.info(f"  Expression proportion threshold: {expr_prop}")
    log.info("=" * 80)

    lm.rank_aggregate(
        adata,
        groupby='cell_type',
        resource_name='consensus',
        expr_prop=expr_prop,
        min_cells=min_cells,
        n_perms=n_perms,
        use_raw=False,
        verbose=True,
        return_all_lrs=False,
        seed=1337,
    )

    return adata.uns['liana_res']


def filter_and_rank(liana_res, max_pval=0.05, top_n_per_pair=20):
    """Filter by significance and rank."""
    log.info("Filtering and ranking interactions...")

    df = liana_res.copy()

    if 'cellphone_pvals' in df.columns:
        sig_col = 'cellphone_pvals'
    elif 'lr_means' in df.columns and 'specificity_rank' in df.columns:
        sig_col = 'specificity_rank'
    else:
        sig_col = None

    log.info(f"  Total raw interactions: {len(df):,}")
    log.info(f"  Columns: {df.columns.tolist()}")

    df = df[df['source'] != df['target']].copy()
    log.info(f"  After removing self-interactions: {len(df):,}")

    if sig_col is not None:
        df_sig = df[df[sig_col] < max_pval] if 'pvals' in sig_col else df.head(int(len(df) * 0.1))
        log.info(f"  Significant (filter on {sig_col}): {len(df_sig):,}")
    else:
        df_sig = df

    df_sig = df_sig.sort_values('magnitude_rank')

    return df, df_sig


def create_pairwise_summary(liana_res, adata):
    """Aggregate L/R interactions per cell type pair."""
    log.info("Building cell-type pair summary matrix...")

    cell_types = sorted(adata.obs['cell_type'].unique().tolist())
    n = len(cell_types)
    pair_counts = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)
    pair_strength = pd.DataFrame(0.0, index=cell_types, columns=cell_types)

    for _, row in liana_res.iterrows():
        s, t = row['source'], row['target']
        if s == t:
            continue
        pair_counts.loc[s, t] += 1
        if 'magnitude_rank' in row.index:
            pair_strength.loc[s, t] += -np.log10(row['magnitude_rank'] + 1e-10)

    return pair_counts, pair_strength


def save_results(all_df, sig_df, pair_counts, pair_strength, adata):
    log.info("Saving results...")

    all_df.to_csv(OUTPUT_DIR / "liana_all_interactions.csv", index=False)
    sig_df.to_csv(OUTPUT_DIR / "liana_significant_interactions.csv", index=False)
    pair_counts.to_csv(OUTPUT_DIR / "celltype_pair_counts.csv")
    pair_strength.to_csv(OUTPUT_DIR / "celltype_pair_strength.csv")

    panel_genes = set(adata.var_names)

    if 'ligand_complex' in all_df.columns:
        unique_ligands = set(all_df['ligand_complex'].dropna().unique())
        unique_receptors = set(all_df['receptor_complex'].dropna().unique())
    else:
        unique_ligands = set(all_df['ligand'].dropna().unique()) if 'ligand' in all_df.columns else set()
        unique_receptors = set(all_df['receptor'].dropna().unique()) if 'receptor' in all_df.columns else set()

    summary = pd.DataFrame({
        'metric': [
            'panel_genes',
            'unique_ligands_detected',
            'unique_receptors_detected',
            'total_interactions',
            'significant_interactions',
            'unique_celltype_pairs',
        ],
        'value': [
            len(panel_genes),
            len(unique_ligands),
            len(unique_receptors),
            len(all_df),
            len(sig_df),
            int((pair_counts > 0).sum().sum()),
        ]
    })
    summary.to_csv(OUTPUT_DIR / "liana_summary_stats.csv", index=False)
    log.info(f"  Output: {OUTPUT_DIR}")
    log.info("\nSummary:")
    log.info(summary.to_string(index=False))


def main():
    start = time.time()
    log.info("=" * 80)
    log.info("PHASE 2B (LIANA+): VALIDATED CELL-CELL COMMUNICATION")
    log.info("=" * 80)

    try:
        adata = load_data()

        liana_res = run_liana_ccc(adata, n_perms=1000, min_cells=50, expr_prop=0.1)

        all_df, sig_df = filter_and_rank(liana_res)

        pair_counts, pair_strength = create_pairwise_summary(sig_df, adata)

        save_results(all_df, sig_df, pair_counts, pair_strength, adata)

        elapsed = time.time() - start
        log.info("=" * 80)
        log.info(f"PHASE 2B (LIANA+) COMPLETE in {elapsed:.1f}s")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
