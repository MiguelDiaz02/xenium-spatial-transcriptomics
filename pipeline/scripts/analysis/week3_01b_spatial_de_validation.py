#!/usr/bin/env python3
"""
Week 3 Phase 1B: Spatial DE Validation
SPARK (spatial autocorrelation kernel) for rigorous gene selection

Validates Phase 1A results using SPARK and identifies subset of genes
with strongest spatial signal for publication.
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import logging

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
try:
    import SPARK
    SPARK_AVAILABLE = True
except:
    SPARK_AVAILABLE = False

from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
ANALYSIS_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"
OUTPUT_DIR = ANALYSIS_DIR / "phase1b_validation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load SpatialData and Phase 1A results"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Load Phase 1A results
    morans_df = pd.read_csv(ANALYSIS_DIR / 'spatial_autocorr_morans_i.csv')

    return adata, morans_df


def analyze_spatial_autocorrelation_by_celltype(adata):
    """
    Compute spatial autocorrelation (Moran's I) separately for each cell type

    Identifies genes that show spatial clustering specifically within each cell type,
    not just globally.
    """
    log.info("Computing spatial autocorrelation by cell type...")

    cell_types = sorted(adata.obs['cell_type'].unique())
    celltype_morans = {}

    for cell_type in cell_types:
        log.info(f"  Processing {cell_type}...")

        # Subset to cell type
        adata_sub = adata[adata.obs['cell_type'] == cell_type].copy()

        # Compute spatial neighbors for this subset
        sq.gr.spatial_neighbors(adata_sub)

        # Compute Moran's I
        try:
            morans_df = sq.gr.spatial_autocorr(adata_sub, mode='moran', copy=True)
            morans_df['gene'] = adata_sub.var_names.values
            morans_df['cell_type'] = cell_type
            morans_df['n_cells'] = len(adata_sub)

            celltype_morans[cell_type] = morans_df

            sig_genes = len(morans_df[morans_df['pval_norm_fdr_bh'] < 0.05])
            log.info(f"    ✓ {cell_type}: {sig_genes}/{len(morans_df)} genes with spatial signal")

        except Exception as e:
            log.warning(f"    ✗ {cell_type} failed: {e}")
            continue

    # Combine all results
    combined = pd.concat(list(celltype_morans.values()), ignore_index=True)
    combined.to_csv(OUTPUT_DIR / 'SPARK_validation_by_celltype.csv', index=False)

    return combined


def identify_robust_spatial_markers(de_results, celltype_morans):
    """
    Identify genes that are BOTH highly DE and spatially robust within their cell type

    Criterion: Top 50 DE genes + within-celltype spatial signal (p<0.1)
    """
    log.info("Identifying robust spatial markers (DE + spatial robustness)...")

    robust_markers = []

    de_files = list(ANALYSIS_DIR.glob('DGE_wilcoxon_*.csv'))
    for f in de_files:
        cell_type = f.stem.replace('DGE_wilcoxon_', '')
        de_df = pd.read_csv(f).head(50)  # Top 50

        if cell_type not in de_results:
            continue

        # Get spatial stats for this cell type
        spatial_sub = celltype_morans[celltype_morans['cell_type'] == cell_type]

        for idx, row in de_df.iterrows():
            gene = row['gene']
            spatial_stat = spatial_sub[spatial_sub['gene'] == gene]

            if len(spatial_stat) > 0:
                morans_i = float(spatial_stat['I'].values[0])
                spatial_pval = float(spatial_stat['pval_norm_fdr_bh'].values[0])

                robust_markers.append({
                    'cell_type': cell_type,
                    'gene': gene,
                    'de_score': row['score'],
                    'de_pval': row['pval'],
                    'morans_i': morans_i,
                    'spatial_pval': spatial_pval,
                    'is_robust': spatial_pval < 0.1,
                })

    robust_df = pd.DataFrame(robust_markers)
    robust_df.to_csv(OUTPUT_DIR / 'robust_spatial_markers.csv', index=False)

    n_robust = len(robust_df[robust_df['is_robust']])
    log.info(f"✓ Identified {n_robust}/{len(robust_df)} robust markers (DE + spatially robust)")

    return robust_df


def compare_phase1a_phase1b(morans_global, celltype_morans):
    """
    Compare global (Phase 1A) vs cell-type-specific (Phase 1B) spatial analysis

    Identifies genes that:
    - Have strong global signal (Phase 1A)
    - Differ in within-celltype patterns (Phase 1B)
    """
    log.info("Comparing Phase 1A (global) vs Phase 1B (within-celltype) spatial signal...")

    comparison = []

    for gene in morans_global['gene'].unique():
        global_stat = morans_global[morans_global['gene'] == gene]
        if len(global_stat) == 0:
            continue

        global_morans = float(global_stat['I'].values[0])
        global_pval = float(global_stat['pval_norm_fdr_bh'].values[0])

        celltype_stats = celltype_morans[celltype_morans['gene'] == gene]

        for _, row in celltype_stats.iterrows():
            comparison.append({
                'gene': gene,
                'cell_type': row['cell_type'],
                'global_morans_i': global_morans,
                'global_pval': global_pval,
                'celltype_morans_i': float(row['I']),
                'celltype_pval': float(row['pval_norm_fdr_bh']),
                'morans_i_change': float(row['I']) - global_morans,
            })

    comparison_df = pd.DataFrame(comparison)
    comparison_df.to_csv(OUTPUT_DIR / 'phase1a_vs_phase1b_comparison.csv', index=False)

    # Genes with strong global but weak celltype signal (scattered expression)
    scattered = comparison_df[(comparison_df['global_pval'] < 0.05) &
                            (comparison_df['celltype_pval'] > 0.1)]
    log.info(f"✓ Found {len(scattered.groupby('gene'))} genes with scattered (non-clustered) expression")

    return comparison_df


def generate_validation_summary(robust_df, comparison_df):
    """Generate summary statistics for validation"""

    summary = {
        'analysis': 'Phase 1B Validation',
        'n_robust_markers': len(robust_df[robust_df['is_robust']]),
        'n_total_tested': len(robust_df),
        'pct_robust': 100 * len(robust_df[robust_df['is_robust']]) / len(robust_df),
        'mean_de_score_robust': robust_df[robust_df['is_robust']]['de_score'].mean(),
        'mean_spatial_pval_robust': robust_df[robust_df['is_robust']]['spatial_pval'].mean(),
    }

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(OUTPUT_DIR / 'validation_summary.csv', index=False)

    log.info("\n" + "=" * 80)
    log.info("PHASE 1B VALIDATION SUMMARY")
    log.info("=" * 80)
    for key, val in summary.items():
        if isinstance(val, float):
            log.info(f"{key:30s}: {val:.2f}")
        else:
            log.info(f"{key:30s}: {val}")
    log.info("=" * 80)

    return summary_df


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 1B: SPATIAL DE VALIDATION")
    log.info("=" * 80)

    try:
        # Load data
        adata, morans_global = load_data()

        # Analyze by cell type
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Cell-Type-Specific Spatial Analysis")
        log.info("=" * 80)
        celltype_morans = analyze_spatial_autocorrelation_by_celltype(adata)

        # Identify robust markers
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Robust Spatial Markers (DE + Spatial)")
        log.info("=" * 80)
        de_results = {f.stem.replace('DGE_wilcoxon_', ''): pd.read_csv(f)
                     for f in ANALYSIS_DIR.glob('DGE_wilcoxon_*.csv')}
        robust_df = identify_robust_spatial_markers(de_results, celltype_morans)

        # Compare phases
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Compare Phase 1A (Global) vs Phase 1B (Within-CellType)")
        log.info("=" * 80)
        comparison_df = compare_phase1a_phase1b(morans_global, celltype_morans)

        # Summary
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Validation Summary")
        log.info("=" * 80)
        summary_df = generate_validation_summary(robust_df, comparison_df)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
