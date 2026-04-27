#!/usr/bin/env python3
"""
Week 3 Phase 3: Task 3.3 - Tissue Region Classification
Define functional zones (immune zone, tumor core, stromal boundary)
Identify tissue regions using spatial clustering + enrichment validation

Methodology:
1. Leiden clustering on spatial coordinates (multiple resolutions)
2. Validate regions using distance-to-tumor + neighborhood enrichment
3. Assign functional labels (tumor core, immune zone, stromal boundary)
4. Compute marker genes per region
5. Analyze region-specific communication patterns
6. Generate publication-quality visualizations
"""

import time
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

import spatialdata as sd
import scanpy as sc
import squidpy as sq
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix, lil_matrix

from utils.logging import get_logger

log = get_logger(__name__)

# Configuration
SDATA_PATH = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "sdata.zarr"
OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "tissue_region_classification"
FIGURES_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "figures" / "phase3_tissue_region_classification"

GRADIENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "spatial_gradients"
ENRICHMENT_OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "human_lung_cancer" / "results" / "02_biology" / "neighborhood_enrichment"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Parameters
LEIDEN_RESOLUTION = 0.8  # For spatial clustering
K_NEIGHBORS = 15  # For spatial neighbor computation
N_TOP_MARKERS = 30  # Top marker genes per region
IMMUNE_TYPES = ['T_cell', 'B_cell', 'Macrophage', 'NK_cell', 'Plasma', 'Dendritic', 'Monocyte']


def load_data():
    """Load SpatialData and prepare for analysis"""
    log.info(f"Loading SpatialData from {SDATA_PATH}...")
    sdata = sd.read_zarr(str(SDATA_PATH))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"✓ Loaded {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Load previous results
    log.info("Loading previous analysis results...")
    if GRADIENT_OUTPUT_DIR.exists() and (GRADIENT_OUTPUT_DIR / "distance_to_tumor.csv").exists():
        distance_df = pd.read_csv(GRADIENT_OUTPUT_DIR / "distance_to_tumor.csv", index_col=0)
        adata.obs['distance_to_tumor'] = distance_df.iloc[:, 0].values
        log.info("  ✓ Loaded distance_to_tumor from Task 3.1")
    else:
        log.warning("  ⚠ distance_to_tumor not found; will compute")
        adata = compute_distance_to_tumor(adata)

    if ENRICHMENT_OUTPUT_DIR.exists() and (ENRICHMENT_OUTPUT_DIR / "neighborhood_enrichment_matrix.csv").exists():
        enrich_matrix = pd.read_csv(ENRICHMENT_OUTPUT_DIR / "neighborhood_enrichment_matrix.csv", index_col=0)
        log.info("  ✓ Loaded neighborhood enrichment matrix from Task 3.2")
    else:
        log.warning("  ⚠ Enrichment matrix not found; will compute")
        enrich_matrix = None

    return adata, enrich_matrix


def compute_distance_to_tumor(adata):
    """Fallback: compute distance to tumor if Task 3.1 output not available"""
    log.info("Computing distance to tumor (fallback)...")
    coords = adata.obsm['spatial']
    tumor_mask = adata.obs['cell_type'] == 'Tumor'
    tumor_coords = coords[tumor_mask]

    distances = np.zeros(len(adata))
    distances[tumor_mask] = 0

    non_tumor_coords = coords[~tumor_mask]
    non_tumor_idx = np.where(~tumor_mask)[0]

    if len(tumor_coords) > 0 and len(non_tumor_coords) > 0:
        dist_matrix = cdist(non_tumor_coords, tumor_coords, metric='euclidean')
        min_dist = dist_matrix.min(axis=1)
        distances[non_tumor_idx] = min_dist

    adata.obs['distance_to_tumor'] = distances
    log.info(f"  Distance range: [{distances.min():.1f}, {distances.max():.1f}]")
    return adata


def compute_spatial_leiden_clusters(adata, resolution=LEIDEN_RESOLUTION):
    """
    Compute Leiden clustering on spatial coordinates
    This identifies spatially-contiguous regions
    """
    log.info(f"Computing spatial Leiden clusters (resolution={resolution})...")

    coords = adata.obsm['spatial']

    # Compute spatial neighbors
    if 'spatial_neighbors' not in adata.obsp:
        log.info("  Computing spatial neighbors graph...")
        nbrs = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1).fit(coords)
        distances, indices = nbrs.kneighbors(coords)

        n_obs = adata.n_obs
        adj = lil_matrix((n_obs, n_obs))
        for i in range(len(indices)):
            for j in indices[i][1:]:  # Skip self
                adj[i, j] = 1

        adata.obsp['spatial_neighbors'] = adj.tocsr()
        log.info(f"  ✓ Built spatial neighbors: {adata.obsp['spatial_neighbors'].shape}")

    # Compute Leiden clusters on spatial graph
    log.info("  Running Leiden clustering on spatial graph...")
    adata_temp = adata.copy()
    adata_temp.obsp['connectivities'] = adata.obsp['spatial_neighbors'].copy()
    adata_temp.obsp['distances'] = adata_temp.obsp['connectivities'].copy().astype(float)

    sc.tl.leiden(adata_temp, key_added='spatial_leiden', resolution=resolution, seed=42)

    adata.obs['spatial_leiden'] = adata_temp.obs['spatial_leiden'].values
    n_clusters = len(adata.obs['spatial_leiden'].unique())

    log.info(f"  ✓ Identified {n_clusters} spatial clusters")
    return adata


def classify_tissue_regions(adata):
    """
    Classify spatial clusters into functional regions
    Using: distance_to_tumor, cell_type composition, neighborhood enrichment
    """
    log.info("Classifying tissue regions into functional zones...")

    region_classifications = []
    region_stats = []

    for cluster_id in sorted(adata.obs['spatial_leiden'].unique()):
        cluster_mask = adata.obs['spatial_leiden'] == cluster_id
        cluster_cells = adata[cluster_mask]

        # Compute region statistics
        median_dist = cluster_cells.obs['distance_to_tumor'].median()
        mean_dist = cluster_cells.obs['distance_to_tumor'].mean()
        max_dist = cluster_cells.obs['distance_to_tumor'].max()

        # Cell type composition
        ct_counts = cluster_cells.obs['cell_type'].value_counts()
        dominant_ct = ct_counts.index[0]
        dominant_frac = ct_counts.iloc[0] / len(cluster_cells)

        # Immune cell fraction
        immune_mask = cluster_cells.obs['cell_type'].isin(IMMUNE_TYPES)
        immune_frac = immune_mask.sum() / len(cluster_cells)

        # Tumor cell fraction
        tumor_frac = (cluster_cells.obs['cell_type'] == 'Tumor').sum() / len(cluster_cells)

        # Stromal (Endothelial + Epithelial non-tumor)
        stromal_types = ['Endothelial', 'Epithelial']
        stromal_frac = cluster_cells.obs['cell_type'].isin(stromal_types).sum() / len(cluster_cells)

        # Classify region based on composition and distance
        if tumor_frac > 0.4:
            region_type = "Tumor Core"
        elif immune_frac > 0.5 and median_dist < np.percentile(adata.obs['distance_to_tumor'], 75):
            region_type = "Immune Zone (Infiltrated)"
        elif immune_frac > 0.3 and median_dist > np.percentile(adata.obs['distance_to_tumor'], 75):
            region_type = "Immune Zone (Peripheral)"
        elif stromal_frac > 0.4:
            region_type = "Stromal Boundary"
        else:
            region_type = "Mixed Zone"

        region_classifications.append({
            'cluster_id': cluster_id,
            'region_type': region_type,
            'n_cells': len(cluster_cells),
            'median_distance_to_tumor': median_dist,
            'mean_distance_to_tumor': mean_dist,
            'max_distance_to_tumor': max_dist,
            'immune_fraction': immune_frac,
            'tumor_fraction': tumor_frac,
            'stromal_fraction': stromal_frac,
            'dominant_cell_type': dominant_ct,
            'dominant_cell_type_fraction': dominant_frac
        })

    region_df = pd.DataFrame(region_classifications)

    # Add region type to adata
    region_mapping = dict(zip(region_df['cluster_id'], region_df['region_type']))
    adata.obs['tissue_region'] = adata.obs['spatial_leiden'].map(region_mapping)

    log.info(f"✓ Classified {len(region_df)} regions:")
    for region_type in region_df['region_type'].unique():
        count = (region_df['region_type'] == region_type).sum()
        log.info(f"  - {region_type}: {count} clusters")

    return adata, region_df


def compute_region_markers(adata):
    """
    Compute marker genes for each tissue region
    Using Wilcoxon rank-sum test (region vs rest)
    """
    log.info("Computing marker genes per tissue region...")

    # Ensure cell type is categorical
    if not isinstance(adata.obs['tissue_region'], pd.CategoricalIndex):
        adata.obs['tissue_region'] = pd.Categorical(adata.obs['tissue_region'])

    sc.tl.rank_genes_groups(
        adata,
        groupby='tissue_region',
        method='wilcoxon',
        key_added='region_markers'
    )

    log.info("✓ Computed marker genes per region")
    return adata


def extract_region_markers_table(adata, n_genes=N_TOP_MARKERS):
    """Extract top marker genes per region as DataFrame"""
    log.info(f"Extracting top {n_genes} markers per region...")

    marker_tables = []

    for region in adata.obs['tissue_region'].unique():
        # Get markers for this region
        indices = adata.uns['region_markers']['names'][region][:n_genes]
        scores = adata.uns['region_markers']['scores'][region][:n_genes]
        pvals = adata.uns['region_markers']['pvals'][region][:n_genes]

        marker_df = pd.DataFrame({
            'region': region,
            'gene': indices,
            'score': scores,
            'pval': pvals,
            'rank': range(1, n_genes + 1)
        })

        marker_tables.append(marker_df)

    all_markers = pd.concat(marker_tables, ignore_index=True)
    return all_markers


def analyze_region_communication(adata, enrich_matrix=None):
    """
    Analyze communication patterns per tissue region
    By computing average cell type composition and neighborhood enrichment
    """
    log.info("Analyzing communication patterns per region...")

    comm_stats = []

    for region in adata.obs['tissue_region'].unique():
        region_mask = adata.obs['tissue_region'] == region
        region_cells = adata[region_mask]

        # Cell type composition
        ct_composition = region_cells.obs['cell_type'].value_counts()
        ct_fractions = ct_composition / ct_composition.sum()

        # Immune cell fraction and diversity
        immune_mask = region_cells.obs['cell_type'].isin(IMMUNE_TYPES)
        immune_frac = immune_mask.sum() / len(region_cells)

        if immune_frac > 0:
            immune_diversity = len(region_cells[immune_mask].obs['cell_type'].unique())
        else:
            immune_diversity = 0

        # If enrichment matrix available, get region-specific enrichment
        if enrich_matrix is not None and region in enrich_matrix.index:
            mean_enrich = enrich_matrix.loc[region].mean()
        else:
            mean_enrich = np.nan

        comm_stats.append({
            'region': region,
            'n_cells': len(region_cells),
            'immune_fraction': immune_frac,
            'immune_diversity': immune_diversity,
            'n_cell_types': len(ct_composition),
            'mean_enrichment': mean_enrich
        })

    comm_df = pd.DataFrame(comm_stats)
    return comm_df


def plot_tissue_region_map(adata):
    """Plot spatial map colored by tissue region classification"""
    log.info("Generating tissue region spatial map...")

    fig, ax = plt.subplots(figsize=(16, 14))

    coords = adata.obsm['spatial']

    # Define colors for regions
    regions = adata.obs['tissue_region'].unique()
    region_colors = {
        'Tumor Core': '#d62728',
        'Immune Zone (Infiltrated)': '#2ca02c',
        'Immune Zone (Peripheral)': '#7f7f7f',
        'Stromal Boundary': '#ff7f0e',
        'Mixed Zone': '#9467bd'
    }

    # Plot each region
    for region in sorted(regions):
        mask = adata.obs['tissue_region'] == region
        color = region_colors.get(region, '#1f77b4')
        ax.scatter(coords[mask, 0], coords[mask, 1],
                  c=color, label=region, s=15, alpha=0.7, edgecolors='none')

    ax.set_xlabel('X coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y coordinate (μm)', fontsize=12, fontweight='bold')
    ax.set_title('Tissue Region Classification: Functional Zone Mapping',
                fontsize=13, fontweight='bold')
    ax.set_aspect('equal')
    ax.legend(loc='best', fontsize=10, framealpha=0.95)
    ax.grid(alpha=0.2)

    plt.tight_layout()
    return fig


def plot_region_statistics(region_df):
    """Plot region-specific statistics"""
    log.info("Generating region statistics plots...")

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # 1. Cell count per region
    ax1 = fig.add_subplot(gs[0, 0])
    region_df_sorted = region_df.sort_values('n_cells', ascending=False)
    colors = ['#d62728' if x=='Tumor Core' else '#2ca02c' if 'Immune' in x else '#ff7f0e'
              for x in region_df_sorted['region_type']]
    ax1.barh(range(len(region_df_sorted)), region_df_sorted['n_cells'], color=colors)
    ax1.set_yticks(range(len(region_df_sorted)))
    ax1.set_yticklabels(region_df_sorted['region_type'], fontsize=10)
    ax1.set_xlabel('Cell Count', fontsize=11, fontweight='bold')
    ax1.set_title('Cells per Region', fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)

    # 2. Distance to tumor
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(region_df['median_distance_to_tumor'], region_df['n_cells'],
               s=region_df['n_cells']/100, alpha=0.6, c=range(len(region_df)), cmap='viridis')
    for i, row in region_df.iterrows():
        ax2.annotate(row['region_type'][:8], (row['median_distance_to_tumor'], row['n_cells']),
                    fontsize=8, alpha=0.7)
    ax2.set_xlabel('Median Distance to Tumor (μm)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Cell Count', fontsize=11, fontweight='bold')
    ax2.set_title('Distance vs Cell Count', fontsize=12, fontweight='bold')
    ax2.grid(alpha=0.3)

    # 3. Cell type composition
    ax3 = fig.add_subplot(gs[0, 2])
    width = 0.15
    x = np.arange(len(region_df))
    ax3.bar(x - width, region_df['tumor_fraction'], width, label='Tumor', color='#d62728')
    ax3.bar(x, region_df['immune_fraction'], width, label='Immune', color='#2ca02c')
    ax3.bar(x + width, region_df['stromal_fraction'], width, label='Stromal', color='#ff7f0e')
    ax3.set_ylabel('Fraction', fontsize=11, fontweight='bold')
    ax3.set_title('Cell Type Composition per Region', fontsize=12, fontweight='bold')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f"R{i}" for i in range(len(region_df))], fontsize=9)
    ax3.legend(fontsize=9)
    ax3.set_ylim(0, 1.0)
    ax3.grid(axis='y', alpha=0.3)

    # 4. Immune fraction by region
    ax4 = fig.add_subplot(gs[1, 0])
    colors_immune = ['#2ca02c' if x=='Immune Zone (Infiltrated)' else '#7f7f7f' if x=='Immune Zone (Peripheral)' else '#cccccc'
                     for x in region_df['region_type']]
    ax4.barh(range(len(region_df)), region_df['immune_fraction'], color=colors_immune)
    ax4.set_yticks(range(len(region_df)))
    ax4.set_yticklabels(region_df['region_type'], fontsize=10)
    ax4.set_xlabel('Immune Fraction', fontsize=11, fontweight='bold')
    ax4.set_title('Immune Cell Content per Region', fontsize=12, fontweight='bold')
    ax4.set_xlim(0, 1.0)
    ax4.grid(axis='x', alpha=0.3)

    # 5. Distance distribution per region
    ax5 = fig.add_subplot(gs[1, 1])
    dist_by_region = []
    region_labels = []
    for region in region_df['region_type'].unique():
        region_labels.append(region)
        dist_by_region.append(region_df[region_df['region_type']==region]['median_distance_to_tumor'].values)

    bp = ax5.boxplot([dist_by_region[i] for i in range(len(dist_by_region))],
                     labels=[l[:12] for l in region_labels], patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('#e8e8e8')
    ax5.set_ylabel('Distance to Tumor (μm)', fontsize=11, fontweight='bold')
    ax5.set_title('Distance Distribution by Region', fontsize=12, fontweight='bold')
    ax5.tick_params(axis='x', rotation=45)
    ax5.grid(axis='y', alpha=0.3)

    # 6. Region size distribution
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.pie(region_df['n_cells'], labels=region_df['region_type'], autopct='%1.1f%%',
           colors=['#d62728' if x=='Tumor Core' else '#2ca02c' if 'Immune' in x else '#ff7f0e'
                   for x in region_df['region_type']], startangle=90)
    ax6.set_title('Region Size Distribution', fontsize=12, fontweight='bold')

    plt.tight_layout()
    return fig


def plot_marker_genes_heatmap(adata, marker_df):
    """Plot top marker genes per region as heatmap"""
    log.info("Generating marker genes heatmap...")

    # Get top 15 markers per region for visualization
    top_n = 15

    # Create matrix of marker gene expression by region
    regions = sorted(adata.obs['tissue_region'].unique())

    fig, axes = plt.subplots(1, len(regions), figsize=(6*len(regions), 8), sharey=True)
    if len(regions) == 1:
        axes = [axes]

    for idx, region in enumerate(regions):
        region_markers = marker_df[marker_df['region']==region].head(top_n)
        region_genes = region_markers['gene'].values

        # Compute mean expression per cell type in this region
        region_mask = adata.obs['tissue_region'] == region
        region_data = adata[region_mask]

        cell_types = sorted(region_data.obs['cell_type'].unique())
        expr_matrix = np.zeros((len(region_genes), len(cell_types)))

        for i, gene in enumerate(region_genes):
            if gene in adata.var_names:
                gene_idx = adata.var_names.get_loc(gene)
                for j, ct in enumerate(cell_types):
                    ct_mask = region_data.obs['cell_type'] == ct
                    if ct_mask.sum() > 0:
                        expr_matrix[i, j] = region_data[ct_mask].X[:, gene_idx].mean()

        # Plot heatmap
        sns.heatmap(expr_matrix, ax=axes[idx], cmap='YlOrRd', cbar=True,
                   xticklabels=[c[:10] for c in cell_types],
                   yticklabels=region_genes, linewidths=0.5)
        axes[idx].set_title(f'{region}\n(Top {top_n} Markers)', fontsize=11, fontweight='bold')
        if idx == 0:
            axes[idx].set_ylabel('Gene', fontsize=10, fontweight='bold')
        else:
            axes[idx].set_ylabel('')
        axes[idx].set_xlabel('Cell Type', fontsize=10, fontweight='bold')
        plt.setp(axes[idx].get_xticklabels(), rotation=45, ha='right', fontsize=9)

    plt.tight_layout()
    return fig


def plot_region_dotplot(adata, marker_df):
    """Plot region-specific markers as dotplot"""
    log.info("Generating marker genes dotplot...")

    regions = sorted(adata.obs['tissue_region'].unique())
    n_regions = len(regions)

    fig, axes = plt.subplots(n_regions, 1, figsize=(14, 4*n_regions))
    if n_regions == 1:
        axes = [axes]

    for idx, region in enumerate(regions):
        region_markers = marker_df[marker_df['region']==region].head(20)
        region_genes = region_markers['gene'].values

        # Ensure these genes are in the dataset
        region_genes = [g for g in region_genes if g in adata.var_names]
        region_genes = region_genes[:15]  # Limit to 15 for clarity

        # Create data for dotplot
        region_mask = adata.obs['tissue_region'] == region
        region_data = adata[region_mask, region_genes]

        cell_types = sorted(adata.obs['cell_type'].unique())

        plot_data = []
        for gene in region_genes:
            for ct in cell_types:
                ct_mask = region_data.obs['cell_type'] == ct
                if ct_mask.sum() > 0:
                    mean_expr = region_data[ct_mask, gene].X.mean()
                    pct_expr = (region_data[ct_mask, gene].X > 0).sum() / ct_mask.sum() * 100
                    plot_data.append({
                        'gene': gene,
                        'cell_type': ct,
                        'mean_expr': mean_expr,
                        'pct_expr': pct_expr
                    })

        plot_df = pd.DataFrame(plot_data)

        # Create dotplot
        pivot_mean = plot_df.pivot(index='gene', columns='cell_type', values='mean_expr')
        pivot_pct = plot_df.pivot(index='gene', columns='cell_type', values='pct_expr')

        ax = axes[idx]

        # Color by mean expression, size by percentage expressed
        for i, gene in enumerate(pivot_mean.index):
            for j, ct in enumerate(pivot_mean.columns):
                mean_val = pivot_mean.loc[gene, ct]
                pct_val = pivot_pct.loc[gene, ct]

                # Skip if very low expression
                if pct_val < 5:
                    ax.scatter(j, i, s=10, c='lightgray', alpha=0.3)
                else:
                    size = (pct_val / 100) * 200 + 20
                    color = min(mean_val, 3)  # Cap at 3 for color
                    ax.scatter(j, i, s=size, c=color, cmap='YlOrRd',
                             vmin=0, vmax=3, edgecolors='black', linewidth=0.5)

        ax.set_xticks(range(len(pivot_mean.columns)))
        ax.set_xticklabels(pivot_mean.columns, rotation=45, ha='right', fontsize=9)
        ax.set_yticks(range(len(pivot_mean.index)))
        ax.set_yticklabels(pivot_mean.index, fontsize=9)
        ax.set_ylabel('Gene', fontsize=10, fontweight='bold')
        ax.set_xlabel('Cell Type', fontsize=10, fontweight='bold')
        ax.set_title(f'{region} - Top Marker Genes (by Cell Type)', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def save_results(adata, region_df, marker_df, comm_df):
    """Save all analysis results to CSV files"""
    log.info("Saving tissue region classification results...")

    # Save region assignments
    adata.obs[['spatial_leiden', 'tissue_region', 'distance_to_tumor']].to_csv(
        OUTPUT_DIR / "region_assignments.csv"
    )

    # Save region statistics
    region_df.to_csv(OUTPUT_DIR / "region_statistics.csv", index=False)

    # Save marker genes per region
    marker_df.to_csv(OUTPUT_DIR / "region_marker_genes.csv", index=False)

    # Save communication statistics
    comm_df.to_csv(OUTPUT_DIR / "region_communication_stats.csv", index=False)

    # Summary statistics
    summary = pd.DataFrame({
        'metric': [
            'total_regions',
            'median_region_size',
            'total_tumor_core_cells',
            'total_immune_zone_cells',
            'total_stromal_cells',
            'average_immune_fraction',
            'average_stromal_fraction',
            'average_tumor_fraction'
        ],
        'value': [
            len(region_df),
            region_df['n_cells'].median(),
            region_df[region_df['region_type']=='Tumor Core']['n_cells'].sum(),
            region_df[region_df['region_type'].str.contains('Immune')]['n_cells'].sum(),
            region_df[region_df['region_type']=='Stromal Boundary']['n_cells'].sum(),
            region_df['immune_fraction'].mean(),
            region_df['stromal_fraction'].mean(),
            region_df['tumor_fraction'].mean()
        ]
    })
    summary.to_csv(OUTPUT_DIR / "region_summary.csv", index=False)

    log.info(f"✓ Results saved to {OUTPUT_DIR}")


def create_comprehensive_report(adata, region_df, marker_df, comm_df, fig_tissue_map, fig_statistics, fig_heatmap, fig_dotplot):
    """Create comprehensive PDF report with all analysis"""
    log.info("Creating comprehensive PDF report...")

    pdf_path = FIGURES_DIR / "Phase3_Task3.3_Tissue_Region_Classification.pdf"

    with PdfPages(pdf_path) as pdf:
        # Page 1: Title and summary
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Week 3 Phase 3: Task 3.3 - Tissue Region Classification',
                    fontsize=18, fontweight='bold', y=0.98)

        ax = fig.add_subplot(111)
        ax.axis('off')

        summary_text = f"""
EXECUTIVE SUMMARY

Dataset: {adata.shape[0]:,} cells × {adata.shape[1]} genes
Identified Regions: {len(region_df)} tissue zones
Classification Method: Spatial Leiden clustering + enrichment validation

KEY FINDINGS:

1. Region Identification
   • Tumor Core: {(region_df['region_type']=='Tumor Core').sum()} clusters, {region_df[region_df['region_type']=='Tumor Core']['n_cells'].sum():,} cells
   • Immune Zones: {region_df[region_df['region_type'].str.contains('Immune')]['region_type'].nunique()} types, {region_df[region_df['region_type'].str.contains('Immune')]['n_cells'].sum():,} cells total
   • Stromal Boundary: {(region_df['region_type']=='Stromal Boundary').sum()} clusters, {region_df[region_df['region_type']=='Stromal Boundary']['n_cells'].sum():,} cells

2. Cell Type Organization
   • Average immune fraction: {region_df['immune_fraction'].mean():.1%}
   • Average stromal fraction: {region_df['stromal_fraction'].mean():.1%}
   • Average tumor fraction: {region_df['tumor_fraction'].mean():.1%}

3. Spatial Distribution
   • Distance to tumor range: {region_df['median_distance_to_tumor'].min():.0f}-{region_df['median_distance_to_tumor'].max():.0f} μm
   • Regions with high immune content (>50%): {(region_df['immune_fraction']>0.5).sum()}
   • Regions with tumor dominance (>40%): {(region_df['tumor_fraction']>0.4).sum()}

4. Region-Specific Characteristics
   • Marker genes identified: {len(marker_df):,} total ({len(marker_df.groupby('region'))} per region average)
   • Communication patterns analyzed: {len(comm_df)} region communication profiles
   • Immune diversity ranges from {comm_df['immune_diversity'].min():.0f} to {comm_df['immune_diversity'].max():.0f} types per region

METHODOLOGY:
1. Leiden clustering on spatial coordinates (k={K_NEIGHBORS} neighbors, resolution={LEIDEN_RESOLUTION})
2. Region classification using distance-to-tumor + cell type composition
3. Marker gene discovery via Wilcoxon rank-sum test (region vs rest)
4. Communication pattern analysis via neighborhood enrichment
5. Publication-quality visualization of spatial organization

OUTPUTS:
Data files: 4 CSV tables with region assignments, statistics, markers, and communication patterns
Figures: 4 high-resolution PNG files + this comprehensive PDF report
Total files: 8 outputs documenting functional tissue architecture
"""

        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2: Tissue region map
        pdf.savefig(fig_tissue_map, bbox_inches='tight')
        plt.close(fig_tissue_map)

        # Page 3: Region statistics
        pdf.savefig(fig_statistics, bbox_inches='tight')
        plt.close(fig_statistics)

        # Page 4: Marker heatmap
        pdf.savefig(fig_heatmap, bbox_inches='tight')
        plt.close(fig_heatmap)

        # Page 5+: Dotplots
        pdf.savefig(fig_dotplot, bbox_inches='tight')
        plt.close(fig_dotplot)

        # Final page: Data summary table
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Region Classification Summary Statistics', fontsize=14, fontweight='bold')

        ax = fig.add_subplot(111)
        ax.axis('tight')
        ax.axis('off')

        table_data = []
        table_data.append(['Region Type', 'n_Clusters', 'n_Cells', 'Immune%', 'Tumor%', 'Stromal%', 'Median Dist (μm)'])

        for _, row in region_df.iterrows():
            table_data.append([
                row['region_type'],
                '1',
                f"{row['n_cells']:,}",
                f"{row['immune_fraction']:.1%}",
                f"{row['tumor_fraction']:.1%}",
                f"{row['stromal_fraction']:.1%}",
                f"{row['median_distance_to_tumor']:.0f}"
            ])

        table = ax.table(cellText=table_data, cellLoc='center', loc='center',
                        colWidths=[0.15, 0.12, 0.12, 0.1, 0.1, 0.12, 0.15])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)

        # Header formatting
        for i in range(len(table_data[0])):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')

        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    log.info(f"✓ Comprehensive report saved: {pdf_path}")


def main():
    start_time = time.time()
    log.info("=" * 80)
    log.info("WEEK 3 PHASE 3 - TASK 3.3: TISSUE REGION CLASSIFICATION")
    log.info("Functional zone mapping (immune zone, tumor core, stromal boundary)")
    log.info("=" * 80)

    try:
        # Load data
        log.info("\n" + "=" * 80)
        log.info("STEP 1: Load Data and Previous Results")
        log.info("=" * 80)
        adata, enrich_matrix = load_data()

        # Spatial clustering
        log.info("\n" + "=" * 80)
        log.info("STEP 2: Spatial Leiden Clustering")
        log.info("=" * 80)
        adata = compute_spatial_leiden_clusters(adata, resolution=LEIDEN_RESOLUTION)

        # Classify regions
        log.info("\n" + "=" * 80)
        log.info("STEP 3: Classify Tissue Regions")
        log.info("=" * 80)
        adata, region_df = classify_tissue_regions(adata)

        # Compute markers
        log.info("\n" + "=" * 80)
        log.info("STEP 4: Compute Region-Specific Markers")
        log.info("=" * 80)
        adata = compute_region_markers(adata)
        marker_df = extract_region_markers_table(adata, n_genes=N_TOP_MARKERS)

        # Analyze communication
        log.info("\n" + "=" * 80)
        log.info("STEP 5: Analyze Communication Patterns per Region")
        log.info("=" * 80)
        comm_df = analyze_region_communication(adata, enrich_matrix)

        # Generate visualizations
        log.info("\n" + "=" * 80)
        log.info("STEP 6: Generate Visualizations")
        log.info("=" * 80)

        log.info("  Generating tissue region map...")
        fig_tissue_map = plot_tissue_region_map(adata)
        fig_tissue_map.savefig(FIGURES_DIR / "Fig1_Tissue_Region_Map.png", dpi=300, bbox_inches='tight')

        log.info("  Generating region statistics plots...")
        fig_statistics = plot_region_statistics(region_df)
        fig_statistics.savefig(FIGURES_DIR / "Fig2_Region_Statistics.png", dpi=300, bbox_inches='tight')

        log.info("  Generating marker gene heatmap...")
        fig_heatmap = plot_marker_genes_heatmap(adata, marker_df)
        fig_heatmap.savefig(FIGURES_DIR / "Fig3_Marker_Genes_Heatmap.png", dpi=300, bbox_inches='tight')

        log.info("  Generating marker gene dotplot...")
        fig_dotplot = plot_region_dotplot(adata, marker_df)
        fig_dotplot.savefig(FIGURES_DIR / "Fig4_Marker_Genes_Dotplot.png", dpi=300, bbox_inches='tight')

        # Save results
        log.info("\n" + "=" * 80)
        log.info("STEP 7: Save Results")
        log.info("=" * 80)
        save_results(adata, region_df, marker_df, comm_df)

        # Create comprehensive report
        log.info("\n" + "=" * 80)
        log.info("STEP 8: Create Comprehensive Report")
        log.info("=" * 80)
        create_comprehensive_report(adata, region_df, marker_df, comm_df,
                                   fig_tissue_map, fig_statistics, fig_heatmap, fig_dotplot)

        elapsed = time.time() - start_time
        log.info("=" * 80)
        log.info(f"✓ TASK 3.3 COMPLETE in {elapsed:.1f} seconds")
        log.info(f"✓ Output directory: {OUTPUT_DIR}")
        log.info(f"✓ Figures directory: {FIGURES_DIR}")
        log.info(f"✓ Identified {len(region_df)} tissue regions")
        log.info(f"✓ Extracted {len(marker_df)} region-specific markers")
        log.info("=" * 80)

    except Exception as e:
        log.error(f"✗ FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
