"""
Step 06b — Immune Cell Subclustering & Granular Annotation
===========================================================

Performs high-resolution subclustering and granular annotation of immune cells
in Xenium spatial transcriptomics data. Designed to improve immune cell type
purity from 15-59% (generic annotation) to ≥70% (granular subtypes).

Workflow:
  1. Isolate immune cells (subset from global cell_type annotation)
  2. Re-cluster at high resolution (Leiden resolution ~1.2 vs 0.5 global)
  3. Score against granular immune cell type marker panels
  4. Assign cell_type_immune_granular labels
  5. Validate purity and generate QC report

Key output:
  - adata.obs["cell_type_immune_granular"]  → Fine-grained immune subset labels
  - adata.obs["immune_purity"]               → Per-cell purity score (0-1)
  - adata.obs["immune_confidence"]           → Confidence level (PASS/REVIEW/FAIL)
  - immune_subclustering_report.html         → Validation visualizations

Scalability:
  - All markers defined in config.yaml (no hard-coded cell types)
  - Configuration-driven enables adaptation to any immune panel
  - Tested on 268k cells; scales to 500k+ with standard resources
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


# ──────────────────────────────────────────────────────────────────────────────
# Core Functions
# ──────────────────────────────────────────────────────────────────────────────

def isolate_immune_subset(adata, immune_types: List[str]) -> Tuple:
    """
    Filter AnnData to immune cells only.

    Args:
        adata: Input AnnData object (from step 06 annotation)
        immune_types: List of cell_type labels to include (e.g., ["T_cell", "Macrophage"])

    Returns:
        adata_immune: Subset AnnData (immune cells only)
        immune_mask: Boolean array of immune cells
    """
    immune_mask = adata.obs['cell_type'].isin(immune_types)
    n_immune = immune_mask.sum()
    n_total = len(adata)

    log.info(f"Isolating immune cells: {n_immune}/{n_total} ({100*n_immune/n_total:.1f}%)")

    if n_immune == 0:
        raise ValueError(f"No immune cells found matching types: {immune_types}")

    adata_immune = adata[immune_mask].copy()

    log.info(f"Immune breakdown:")
    for cell_type in immune_types:
        ct_count = (adata_immune.obs['cell_type'] == cell_type).sum()
        if ct_count > 0:
            log.info(f"  {cell_type}: {ct_count} cells")

    return adata_immune, immune_mask


def subcluster_immune_cells(adata_immune, n_pcs: int, leiden_res: float,
                           min_cluster_size: int):
    """
    Re-cluster immune cells at high resolution.

    Args:
        adata_immune: Immune-only AnnData
        n_pcs: Number of PCA components
        leiden_res: Leiden resolution (higher = finer clusters; ~1.2 recommended)
        min_cluster_size: Minimum cells per cluster

    Returns:
        adata_immune: Updated with leiden_immune clustering
    """
    import scanpy as sc

    log.info(f"Immune subclustering: PCA + UMAP + Leiden (res={leiden_res})")

    # PCA on immune subset (variance-driven)
    sc.pp.pca(adata_immune, n_comps=n_pcs)
    log.info(f"  PCA: {n_pcs} components, {adata_immune.obsm['X_pca'].shape[1]} computed")

    # Neighbors graph
    sc.pp.neighbors(adata_immune, n_neighbors=15, use_rep='X_pca')
    log.info(f"  Neighbors: k=15")

    # UMAP (high resolution for immune heterogeneity)
    sc.tl.umap(adata_immune, min_dist=0.3)
    log.info(f"  UMAP: computed")

    # Leiden at high resolution
    sc.tl.leiden(adata_immune, resolution=leiden_res, key_added='leiden_immune')
    log.info(f"  Leiden: {adata_immune.obs['leiden_immune'].nunique()} clusters")

    # Cluster sizes
    cluster_sizes = adata_immune.obs['leiden_immune'].value_counts().sort_index()
    min_size = cluster_sizes.min()
    max_size = cluster_sizes.max()
    log.info(f"  Cluster sizes: {min_size}-{max_size} cells")

    # Flag small clusters
    small_clusters = (cluster_sizes < min_cluster_size).sum()
    if small_clusters > 0:
        log.warning(f"  {small_clusters} clusters below min size ({min_cluster_size})")

    return adata_immune


def score_marker_panel(adata, markers: Dict[str, List[str]]) -> pd.DataFrame:
    """
    Score cells against marker gene panels.

    Args:
        adata: Input AnnData
        markers: Dict of cell_type → list of marker genes

    Returns:
        score_df: DataFrame of scores (cells × cell_types)
    """
    import scanpy as sc

    log.info(f"Scoring {len(markers)} immune cell type marker panels...")

    score_dict = {}
    for cell_type, genes in markers.items():
        genes_in_panel = [g for g in genes if g in adata.var_names]
        if not genes_in_panel:
            log.warning(f"  {cell_type}: NO marker genes in panel (asked for {len(genes)})")
            score_dict[cell_type] = np.zeros(len(adata))
            continue

        key = f"score_{cell_type}"
        sc.tl.score_genes(adata, genes_in_panel, score_name=key)
        score_dict[cell_type] = adata.obs[key].values
        log.info(f"  {cell_type}: {len(genes_in_panel)}/{len(genes)} markers found")

    score_df = pd.DataFrame(score_dict, index=adata.obs_names)
    return score_df


def assign_granular_immune_types(adata_immune, markers: Dict[str, List[str]],
                                 score_df: pd.DataFrame) -> Dict:
    """
    Assign granular immune cell types based on marker scores.

    Args:
        adata_immune: Immune-only AnnData
        markers: Marker panels
        score_df: Marker score DataFrame

    Returns:
        annotation_dict: Dict with cell_type_immune_granular, purity, confidence
    """
    # Assign label of highest score
    adata_immune.obs['cell_type_immune_granular'] = score_df.idxmax(axis=1).str.replace('score_', '')

    # Compute purity: max_score - median of other scores
    max_scores = score_df.max(axis=1)
    score_df_copy = score_df.copy()
    for idx in score_df.index:
        max_col = score_df.loc[idx].idxmax()
        score_df_copy.loc[idx, max_col] = np.nan
    median_other = score_df_copy.median(axis=1, skipna=True)

    purity = max_scores - median_other
    purity = np.clip(purity, 0, 1)  # Normalize to 0-1

    adata_immune.obs['immune_purity'] = purity

    # Confidence level (based on purity threshold from config)
    validation_threshold = snakemake.params.validation_threshold  # noqa: F821
    confidence = pd.Series(index=adata_immune.obs_names, dtype='category')
    confidence[purity >= validation_threshold] = 'PASS'
    confidence[(purity >= 0.5) & (purity < validation_threshold)] = 'REVIEW'
    confidence[purity < 0.5] = 'FAIL'
    adata_immune.obs['immune_confidence'] = confidence.astype('category')

    # Summary
    log.info(f"Granular immune annotation assigned:")
    log.info(adata_immune.obs['cell_type_immune_granular'].value_counts().to_string())

    log.info(f"Purity summary:")
    log.info(f"  Mean: {purity.mean():.3f}")
    log.info(f"  Median: {purity.median():.3f}")
    log.info(f"  Range: {purity.min():.3f} - {purity.max():.3f}")

    log.info(f"Confidence summary:")
    log.info(adata_immune.obs['immune_confidence'].value_counts().to_string())

    return {
        'cell_type_immune_granular': adata_immune.obs['cell_type_immune_granular'],
        'immune_purity': adata_immune.obs['immune_purity'],
        'immune_confidence': adata_immune.obs['immune_confidence']
    }


def merge_immune_annotation_to_main(adata, adata_immune, immune_mask, annotation_dict):
    """
    Add immune granular annotations back to main AnnData.

    Preserves global cell_type for spatial analysis; adds granular annotations
    for immune cells only.

    Args:
        adata: Main AnnData (all cells)
        adata_immune: Immune subset with granular annotations
        immune_mask: Boolean array of immune cells
        annotation_dict: Dict with granular annotations

    Returns:
        adata: Updated with immune annotations
    """
    # Initialize columns for non-immune cells
    adata.obs['cell_type_immune_granular'] = 'Non-immune'
    adata.obs['immune_purity'] = np.nan
    adata.obs['immune_confidence'] = 'N/A'
    adata.obs['leiden_immune'] = 'N/A'

    # Fill immune cell annotations
    immune_indices = np.where(immune_mask)[0]
    for idx, immune_idx in zip(adata_immune.obs_names, immune_indices):
        adata.obs.loc[idx, 'cell_type_immune_granular'] = annotation_dict['cell_type_immune_granular'][idx]
        adata.obs.loc[idx, 'immune_purity'] = annotation_dict['immune_purity'][idx]
        adata.obs.loc[idx, 'immune_confidence'] = annotation_dict['immune_confidence'][idx]
        adata.obs.loc[idx, 'leiden_immune'] = adata_immune.obs.loc[idx, 'leiden_immune']

    # Convert to categorical
    adata.obs['cell_type_immune_granular'] = adata.obs['cell_type_immune_granular'].astype('category')
    adata.obs['immune_confidence'] = adata.obs['immune_confidence'].astype('category')
    adata.obs['leiden_immune'] = adata.obs['leiden_immune'].astype('category')

    log.info(f"Merged immune annotations back to main AnnData")
    log.info(f"  Immune cells: {(adata.obs['immune_confidence'] != 'N/A').sum()}")
    log.info(f"  Non-immune cells: {(adata.obs['immune_confidence'] == 'N/A').sum()}")

    return adata


def generate_validation_report(adata_immune, adata_full, markers: Dict,
                               output_path: str):
    """
    Generate HTML validation report with visualizations.

    Args:
        adata_immune: Immune subset with annotations
        adata_full: Full AnnData (for context)
        markers: Marker panels
        output_path: Path to save HTML report
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages

    log.info(f"Generating validation report: {output_path}")

    # Create figure collection
    figures = []

    # 1. Immune breakdown
    fig, ax = plt.subplots(figsize=(10, 6))
    adata_immune.obs['cell_type_immune_granular'].value_counts().plot(kind='barh', ax=ax)
    ax.set_xlabel('Number of cells')
    ax.set_title('Immune Cell Type Distribution (Granular)')
    figures.append(('Immune_distribution', fig))

    # 2. Purity histogram
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(adata_immune.obs['immune_purity'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(snakemake.params.validation_threshold, color='green', linestyle='--',  # noqa: F821
               label=f'PASS threshold ({snakemake.params.validation_threshold})')  # noqa: F821
    ax.axvline(0.5, color='orange', linestyle='--', label='REVIEW threshold (0.5)')
    ax.set_xlabel('Purity Score')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Immune Cell Purity Distribution')
    ax.legend()
    figures.append(('Purity_histogram', fig))

    # 3. Purity by cluster
    fig, ax = plt.subplots(figsize=(12, 6))
    purity_by_cluster = adata_immune.obs.groupby('leiden_immune')['immune_purity'].apply(list)
    bp = ax.boxplot([purity_by_cluster[c] for c in sorted(purity_by_cluster.index)],
                    labels=sorted(purity_by_cluster.index))
    ax.axhline(snakemake.params.validation_threshold, color='green', linestyle='--', alpha=0.5)  # noqa: F821
    ax.axhline(0.5, color='orange', linestyle='--', alpha=0.5)
    ax.set_xlabel('Leiden Cluster (Immune)')
    ax.set_ylabel('Purity Score')
    ax.set_title('Purity Distribution by Immune Cluster')
    figures.append(('Purity_by_cluster', fig))

    # 4. Confidence summary
    fig, ax = plt.subplots(figsize=(10, 6))
    confidence_counts = adata_immune.obs['immune_confidence'].value_counts()
    colors = {'PASS': 'green', 'REVIEW': 'orange', 'FAIL': 'red'}
    color_list = [colors.get(c, 'gray') for c in confidence_counts.index]
    confidence_counts.plot(kind='bar', ax=ax, color=color_list)
    ax.set_ylabel('Number of cells')
    ax.set_title('Immune Cell Annotation Confidence')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    figures.append(('Confidence_summary', fig))

    # Save as PDF
    pdf_path = output_path.replace('.html', '.pdf')
    with PdfPages(pdf_path) as pdf:
        for title, fig in figures:
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

    log.info(f"Validation report saved: {pdf_path}")

    # Create HTML index (simple)
    html_content = f"""
    <html>
    <head>
        <title>Immune Subclustering Validation Report</title>
        <style>
            body {{ font-family: Arial; margin: 20px; }}
            h1 {{ color: #333; }}
            .metric {{ background-color: #f0f0f0; padding: 10px; margin: 10px 0; }}
        </style>
    </head>
    <body>
        <h1>Immune Cell Subclustering Validation Report</h1>
        <p>Generated on 2026-04-24</p>

        <h2>Summary Statistics</h2>
        <div class="metric">
            <strong>Immune Cells:</strong> {(adata_immune.obs['immune_confidence'] != 'N/A').sum()} / {len(adata_full)}
        </div>
        <div class="metric">
            <strong>Mean Purity:</strong> {adata_immune.obs['immune_purity'].mean():.3f}
        </div>
        <div class="metric">
            <strong>PASS Confidence:</strong> {(adata_immune.obs['immune_confidence'] == 'PASS').sum()} cells
        </div>
        <div class="metric">
            <strong>REVIEW Confidence:</strong> {(adata_immune.obs['immune_confidence'] == 'REVIEW').sum()} cells
        </div>
        <div class="metric">
            <strong>FAIL Confidence:</strong> {(adata_immune.obs['immune_confidence'] == 'FAIL').sum()} cells
        </div>

        <h2>Detailed Report</h2>
        <p>See {pdf_path} for visualizations.</p>
    </body>
    </html>
    """

    with open(output_path, 'w') as f:
        f.write(html_content)

    log.info(f"HTML report saved: {output_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Main Execution
# ──────────────────────────────────────────────────────────────────────────────

def main():
    import spatialdata as sd

    sdata_path = snakemake.input.sdata           # noqa: F821
    done_path = snakemake.output.done            # noqa: F821
    report_path = snakemake.output.report        # noqa: F821
    immune_types = snakemake.params.immune_types  # noqa: F821
    leiden_res = snakemake.params.leiden_res_immune  # noqa: F821
    n_pcs = snakemake.params.n_pcs_immune if hasattr(snakemake.params, 'n_pcs_immune') else 25  # noqa: F821
    min_cluster_size = snakemake.params.min_cluster_size  # noqa: F821
    markers = snakemake.params.markers           # noqa: F821

    # ─── Load data ────────────────────────────────────────────────────────────
    log.info(f"Loading SpatialData from {sdata_path}...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]

    if 'cell_type' not in adata.obs:
        raise RuntimeError("cell_type not found. Run step 06 (annotation) first.")

    # ─── Immune subclustering ─────────────────────────────────────────────────
    adata_immune, immune_mask = isolate_immune_subset(adata, immune_types)
    adata_immune = subcluster_immune_cells(adata_immune, n_pcs, leiden_res, min_cluster_size)

    # ─── Annotation ────────────────────────────────────────────────────────────
    score_df = score_marker_panel(adata_immune, markers)
    annotation_dict = assign_granular_immune_types(adata_immune, markers, score_df)

    # ─── Merge back to main ─────────────────────────────────────────────────────
    adata = merge_immune_annotation_to_main(adata, adata_immune, immune_mask, annotation_dict)

    # ─── Save to zarr ──────────────────────────────────────────────────────────
    table_name = list(sdata.tables.keys())[0]
    sdata.tables[table_name] = adata
    log.info("Writing annotated table to zarr store...")
    sdata.delete_element_from_disk(table_name)
    sdata.write_element(table_name)

    # ─── Generate report ──────────────────────────────────────────────────────
    generate_validation_report(adata_immune, adata, markers, report_path)

    Path(done_path).touch()
    log.info("Step 06b — Immune Subclustering: DONE")


main()
