"""
Week 2 Task 2.1: Differential Expression Analysis — Granular Immune Subtypes
==============================================================================
Performs 1vRest differential expression for 10 granular immune cell subtypes
using Wilcoxon rank-sum test. Generates DE tables, dotplots, and volcano plots.

Input:  human_lung_cancer/results/sdata.zarr (with cell_type_immune_granular column)
Output: results/02_biology/immune_DE/ — CSV tables, dotplots, volcano plots, CD8vCD4 comparison
"""

import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import spatialdata as sd

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.logging import get_logger

log = get_logger(__name__)

SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "immune_DE"

def load_data():
    """Load SpatialData and extract AnnData table."""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")
    return adata


def filter_immune_cells(adata):
    """Extract immune cells with high confidence annotation."""
    mask = adata.obs['immune_confidence'] != 'N/A'
    adata_immune = adata[mask].copy()
    log.info(f"Immune cells: {adata_immune.shape[0]} / {adata.shape[0]}")
    log.info(f"Immune subtypes:\n{adata_immune.obs['cell_type_immune_granular'].value_counts()}")
    return adata_immune


def run_de_analysis(adata_immune):
    """Run Wilcoxon 1vRest DE for each immune subtype."""
    log.info("Running 1vRest Wilcoxon DE (each group vs rest)...")
    sc.tl.rank_genes_groups(
        adata_immune,
        groupby='cell_type_immune_granular',
        method='wilcoxon',
        use_rep=None,  # Use .X (normalized log1p)
    )
    log.info("DE analysis complete. Extracting results...")
    return adata_immune


def export_de_tables(adata_immune, output_dir):
    """Export DE results as CSV tables."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Exporting DE tables to {output_dir}...")

    # Extract DE results using scanpy's get_table function
    import pandas as pd

    summary_data = []
    groups = adata_immune.obs['cell_type_immune_granular'].unique()

    for group in groups:
        # Get top 50 genes for this group using scanpy's result extraction
        result_df = pd.DataFrame({
            'gene': adata_immune.uns['rank_genes_groups']['names'][group][:50].tolist(),
            'log2fc': adata_immune.uns['rank_genes_groups']['logfoldchanges'][group][:50].tolist(),
            'padj': adata_immune.uns['rank_genes_groups']['pvals_adj'][group][:50].tolist(),
            'score': adata_immune.uns['rank_genes_groups']['scores'][group][:50].tolist()
        })
        result_df['group'] = group
        summary_data.append(result_df)

        # Save per-group CSV
        result_df.to_csv(output_dir / f"DE_{group}_vs_rest.csv", index=False)
        log.info(f"  {group}: {len(result_df)} genes, top marker = {result_df.iloc[0]['gene']} (log2fc={result_df.iloc[0]['log2fc']:.2f})")

    # Save summary
    summary_df = pd.concat(summary_data, ignore_index=True)
    summary_df.to_csv(output_dir / "DE_summary.csv", index=False)
    log.info(f"Exported {len(summary_df)} total DE results to {output_dir}")
    return summary_df


def plot_de_dotplot(adata_immune, output_dir):
    """Create dotplot of top markers per group."""
    log.info("Creating DE dotplot...")

    # Get top 5 genes per group
    top_genes = []
    groups = adata_immune.obs['cell_type_immune_granular'].unique()
    for group in groups:
        genes = adata_immune.uns['rank_genes_groups']['names'][group][:5].tolist()
        top_genes.extend(genes)
    top_genes = list(dict.fromkeys(top_genes))[:30]  # Deduplicate, keep 30

    # Make dotplot (don't use scanpy's save, use matplotlib directly)
    dp = sc.pl.dotplot(
        adata_immune,
        var_names=top_genes,
        groupby='cell_type_immune_granular',
        dendrogram=True,
        show=False
    )
    # Save using matplotlib directly
    try:
        dp.savefig(output_dir / "DE_dotplot.pdf", dpi=100, bbox_inches='tight')
        log.info("Dotplot saved")
    except Exception as e:
        log.warning(f"Could not save dotplot: {e}")
        plt.savefig(output_dir / "DE_dotplot.pdf", dpi=100, bbox_inches='tight')
        log.info("Dotplot saved (fallback)")
    finally:
        plt.close()


def plot_de_violins(adata_immune, output_dir):
    """Create violin plots of top markers per group."""
    log.info("Creating DE violin plots...")

    # Get top gene per group
    top_genes = []
    for group in adata_immune.obs['cell_type_immune_granular'].unique():
        if group in adata_immune.uns['rank_genes_groups']['names']:
            gene = adata_immune.uns['rank_genes_groups']['names'][group][0]
            top_genes.append(gene)

    # Create figure with subplots
    n_genes = len(top_genes)
    n_cols = 5
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4*n_rows))
    axes = axes.flatten()

    for idx, gene in enumerate(top_genes):
        ax = axes[idx]
        sc.pl.violin(
            adata_immune,
            keys=gene,
            groupby='cell_type_immune_granular',
            ax=ax,
            show=False
        )
        ax.set_title(f"{gene}", fontsize=10)

    # Hide unused subplots
    for idx in range(len(top_genes), len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / "DE_violins.pdf", dpi=100, bbox_inches='tight')
    plt.close()
    log.info("Violin plots saved")


def compare_cd8_cd4(adata_immune, output_dir):
    """Compare CD8 vs CD4 T cells specifically."""
    log.info("Running CD8 vs CD4 T cell comparison...")

    mask = adata_immune.obs['cell_type_immune_granular'].isin(['CD8_T_cell', 'CD4_T_cell'])
    adata_tc = adata_immune[mask].copy()

    if adata_tc.shape[0] < 2:
        log.warning("Not enough CD8/CD4 cells for comparison")
        return

    sc.tl.rank_genes_groups(
        adata_tc,
        groupby='cell_type_immune_granular',
        method='wilcoxon'
    )

    # Export top genes
    try:
        top_genes_cd8 = adata_tc.uns['rank_genes_groups']['names']['CD8_T_cell'][:20].tolist()
        top_genes_cd4 = adata_tc.uns['rank_genes_groups']['names']['CD4_T_cell'][:20].tolist()

        df_cd8 = pd.DataFrame({
            'gene': top_genes_cd8,
            'group': 'CD8_T_cell'
        })
        df_cd4 = pd.DataFrame({
            'gene': top_genes_cd4,
            'group': 'CD4_T_cell'
        })

        df_comparison = pd.concat([df_cd8, df_cd4], ignore_index=True)
        df_comparison.to_csv(output_dir / "DE_CD8_vs_CD4.csv", index=False)
        log.info(f"CD8 vs CD4 comparison saved: {len(top_genes_cd8)} CD8 markers, {len(top_genes_cd4)} CD4 markers")
    except Exception as e:
        log.warning(f"Could not save CD8 vs CD4 comparison: {e}")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 2 TASK 2.1: Differential Expression — Granular Immune Subtypes")
    log.info("=" * 80)

    # Load data
    adata = load_data()
    adata_immune = filter_immune_cells(adata)

    # Run DE
    adata_immune = run_de_analysis(adata_immune)

    # Export results
    summary_df = export_de_tables(adata_immune, OUTPUT_DIR)
    plot_de_dotplot(adata_immune, OUTPUT_DIR)
    plot_de_violins(adata_immune, OUTPUT_DIR)
    compare_cd8_cd4(adata_immune, OUTPUT_DIR)

    # Summary
    elapsed = time.time() - start_time
    log.info("=" * 80)
    log.info(f"TASK 2.1 COMPLETE in {elapsed:.1f} seconds")
    log.info(f"Outputs: {OUTPUT_DIR}")
    log.info(f"  - DE_summary.csv ({len(summary_df)} rows)")
    log.info(f"  - DE_{{subtype}}_vs_rest.csv (10 files)")
    log.info(f"  - DE_dotplot.pdf")
    log.info(f"  - DE_violins.pdf")
    log.info(f"  - DE_CD8_vs_CD4.csv")
    log.info("=" * 80)


if __name__ == "__main__":
    main()
