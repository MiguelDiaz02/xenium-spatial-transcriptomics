# Immune Cell Annotation Strategy for Xenium Spatial Transcriptomics

**Document Version:** 1.0  
**Date:** 2026-04-24  
**Status:** IMPLEMENTATION PLAN  
**Target Completion:** Week 2 (May 1-5, 2026)

---

## Executive Summary

This document outlines a **rigorous, reproducible workflow for comprehensive immune cell subclustering and annotation** in Xenium spatial transcriptomics data. The strategy addresses the current limitation (immune cell purity 15-59%) by implementing:

1. **Immune cell subset isolation** from the main AnnData object
2. **High-resolution subclustering** (Leiden clustering within immune population only)
3. **Granular immune cell type annotation** using panel-specific markers
4. **Quantitative purity validation** with scoring and heatmap visualization
5. **Snakemake integration** as a new step (06b) for reproducibility and future scalability

This workflow is **critical for future projects** using immunology-focused Xenium panels (dedicated immune gene panels for lung and liver), where immune cell heterogeneity is the primary research question.

---

## Part 1: Scientific Strategy

### 1.1 Problem Statement

**Current State:**
- Global Leiden clustering (15 clusters) produces generic "immune" categories
- Manual marker scoring uses only pan-immune markers (CD45, CD3D, CD68, etc.)
- No resolution of immune cell subtypes: T-cell subsets (CD4+, CD8+, Treg), macrophage polarization (M1/M2), NK maturity states
- Result: **Low purity (15-59%)** because heterogeneous cells are grouped together

**Why This Matters:**
1. **For current lung cancer dataset:** DE analysis, L/R interactions contaminated; spatial patterns masked
2. **For future immune panels:** Fundamental research question is immune heterogeneity → current approach fails
3. **For publication:** Reviewers demand validation of cell type identities, especially immune cells

### 1.2 Solution: Immune Subclustering Workflow

```
Input: Global AnnData (after step 06 annotation)
  ↓
Step 1: Isolate immune cells
  - Filter: cell_type ∈ {T_cell, B_cell, Macrophage, NK_cell, Mast_cell, ...}
  - Subset AnnData → adata_immune
  ↓
Step 2: High-resolution reclustering
  - Recompute PCA on immune subset (variance-driven)
  - Recompute UMAP (better resolution for immune heterogeneity)
  - Leiden clustering with HIGHER resolution (e.g., 1.0-1.5 vs 0.5 global)
  ↓
Step 3: Granular annotation
  - Define immune subset markers (CD4, CD8, FOXP3 for T-cells; etc.)
  - Score against marker panels
  - Assign fine-grained labels: "CD8_T_cell", "CD4_T_cell", "Treg", "M1_Macrophage", etc.
  ↓
Step 4: Purity validation
  - Compute marker consistency heatmap
  - Calculate purity scores per subtype
  - Flag low-confidence clusters for review
  ↓
Step 5: Merge & export
  - Add `obs["cell_type_immune_granular"]` to main AnnData
  - Preserve global annotation for spatial analysis
  - Export QC report with purity metrics
  ↓
Output: Enhanced AnnData with granular immune annotation
```

### 1.3 Immune Cell Markers by Subtype

Based on 289-gene Xenium v1 panel and literature consensus:

#### T Cells & NK Cells
```yaml
CD8_T_cell:
  markers: [CD8A, CD8B, GZMA, GZMB, PRF1]
  validation: [CD3D, CD3E]  # Must express
  exclusion: [CD4, FOXP3]   # Must NOT express
  
CD4_T_cell:
  markers: [CD4, IL7R, TCF7, CCR7]
  validation: [CD3D, CD3E]
  exclusion: [CD8A, FOXP3]
  
Treg:
  markers: [FOXP3, IL2RA, CTLA4, ENTPD1]
  validation: [CD4, CD3D, CD3E]
  exclusion: [CD8A, CD25_low]  # FOXP3+CD4+ only
  
NK_cell:
  markers: [GNLY, NKG7, FGFBP2, FCGR3A]
  validation: [NCAM1]
  exclusion: [CD3D, CD3E, CD4]
```

#### Myeloid Cells
```yaml
M1_Macrophage:
  markers: [IL1B, TNF, IL6, NOS2, CXCL9, CXCL10]
  validation: [CD68, CD163]
  polarization: "pro-inflammatory"
  
M2_Macrophage:
  markers: [IL10, TGM2, MRC1, ARG1, CD163, MSR1]
  validation: [CD68]
  polarization: "anti-inflammatory"
  
Macrophage_tissue:
  markers: [CD68, MARCO, SPI1, MERTK]
  validation: [ITGAM]  # CD11b
  
Monocyte:
  markers: [CD14, ITGAM, LYZ, S100A12]
  exclusion: [CD3D, CD19]
  
Dendritic_cell:
  markers: [CD11C, HLA-DRA, FCER1A, LAMP3]
  validation: [CD1C, CLEC9A]
```

#### B Cells & Plasma Cells
```yaml
B_cell_naive:
  markers: [MS4A1, CD19, IGHM, IGHD]
  exclusion: [IGHG1, IGHG3, CD27_high]
  
B_cell_memory:
  markers: [MS4A1, CD19, IGHG1, IGHG3, CD27]
  exclusion: [IGHM, IGHD]
  
Plasma_cell:
  markers: [IGHG1, IGHG3, XBP1, IRF4, SDC1]
  exclusion: [MS4A1, CD19]
```

**Note:** Markers selected from 289-gene panel. For 5K immune panels, expand with:
- T-cell exhaustion: PDCD1, CTLA4, HAVCR2, TIGIT, LAG3
- T-cell activation: IFNG, IL2, GZMA, GZMB
- Macrophage metabolism: LDHA, PFKFB3 (M1), IDO1, IL4RA (M2)

### 1.4 Validation Metrics

**Purity Score (per cell):**
```
purity = max(score_markers) - median(score_other_types)
```

**Cluster-level validation:**
```
- Median purity ≥ 0.7 → PASS (confident subtype)
- Median purity 0.5-0.7 → REVIEW (ambiguous, flag for inspection)
- Median purity < 0.5 → FAIL (mixed population, may need remerge)
```

**Visualization:**
1. **Heatmap:** Top 5 markers per subtype × all cells, ordered by cluster
2. **QC scatter:** Cell purity (y-axis) vs Leiden cluster (x-axis)
3. **UMAP colored by purity:** Visual detection of low-confidence regions

---

## Part 2: Technical Implementation

### 2.1 New Snakemake Step: 06b_immune_subclustering

**Location in DAG:**
```
Step 05 (reduction) → Step 06 (annotation) → Step 06b (immune_subclustering) → Step 07 (denoising)
```

**Dependencies:**
- Input: `results/sdata.zarr` (after step 06, contains `cell_type` annotation)
- Output: Same `sdata.zarr` with new `obs["cell_type_immune_granular"]` + QC report

**Snakemake Rule (new file: `rules/06b_immune_subclustering.smk`):**

```snakemake
rule immune_subclustering:
    input:
        sdata = f"{SDATA}",
        done = f"{OUTDIR}/06_annotation.done"
    output:
        done = f"{OUTDIR}/06b_immune_subclustering.done",
        report = f"{OUTDIR}/06b_immune_subclustering_report.html"
    log:
        "logs/06b_immune_subclustering.log"
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 300
    conda:
        "envs/xenium_pipeline.yaml"
    params:
        immune_types = config.get("immune_annotation", {}).get("immune_cell_types", []),
        leiden_res_immune = config.get("immune_annotation", {}).get("leiden_resolution_immune", 1.2),
        min_cluster_size = config.get("immune_annotation", {}).get("min_immune_cluster_size", 50),
        markers = config.get("immune_annotation", {}).get("immune_markers", {})
    script:
        "scripts/06b_immune_subclustering.py"
```

**Config additions** (to `config_lung.yaml` and `config_liver.yaml`):

```yaml
# Immune Cell Subclustering (new section)
immune_annotation:
  enabled: true                          # Enable immune subclustering
  immune_cell_types:                     # Cell types to subset for immune analysis
    - "T_cell"
    - "B_cell"
    - "Macrophage"
    - "NK_cell"
    - "Mast_cell"
  leiden_resolution_immune: 1.2          # Higher than global (0.5) for immune heterogeneity
  min_immune_cluster_size: 50            # Minimum cells per cluster
  n_pcs_immune: 25                       # PCA components for immune subset
  immune_markers:                        # Immune subtype markers
    CD8_T_cell: ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1"]
    CD4_T_cell: ["CD4", "IL7R", "TCF7", "CCR7"]
    Treg: ["FOXP3", "IL2RA", "CTLA4"]
    NK_cell: ["GNLY", "NKG7", "FGFBP2"]
    M1_Macrophage: ["IL1B", "TNF", "IL6", "NOS2"]
    M2_Macrophage: ["IL10", "TGM2", "MRC1", "ARG1"]
    B_cell: ["MS4A1", "CD19", "IGHM"]
    Plasma_cell: ["IGHG1", "IGHG3", "XBP1"]
    Monocyte: ["CD14", "LYZ", "ITGAM"]
    Dendritic_cell: ["CD1C", "LAMP3"]
  validation_threshold: 0.7              # Min purity for "PASS" status
```

### 2.2 New Script: `06b_immune_subclustering.py`

**Key functions:**

```python
def isolate_immune_subset(adata, immune_types):
    """Filter to immune cells only."""
    immune_mask = adata.obs['cell_type'].isin(immune_types)
    return adata[immune_mask].copy()

def subcluster_immune_cells(adata_immune, n_pcs, leiden_res, min_cluster_size):
    """Re-cluster immune cells at high resolution."""
    # Recompute PCA (variance-driven for immune subset)
    # Recompute UMAP (better resolution)
    # Leiden with high resolution
    # Return: adata_immune with new 'leiden_immune' cluster assignment

def annotate_immune_subtypes(adata_immune, markers, validation_threshold):
    """Assign granular immune subtypes based on marker scoring."""
    # Score against each marker panel
    # Assign label of highest score
    # Compute purity metrics
    # Flag low-confidence cells
    # Return: adata_immune with 'cell_type_immune_granular' + 'purity' columns

def validate_immune_purity(adata_immune, markers):
    """Generate validation heatmap and QC metrics."""
    # Marker consistency heatmap (top markers × cells, ordered by cluster)
    # Purity histogram
    # Per-cluster purity boxplot
    # Return: HTML report

def merge_immune_annotation_to_main(adata, adata_immune):
    """Add immune granular annotation back to main AnnData."""
    # Merge cell_type_immune_granular column
    # Preserve global cell_type for spatial analysis
    # Return: adata with both annotations
```

---

## Part 3: Execution Plan (Week 2)

### Timeline

| Date | Task | Hours | Output |
|------|------|-------|--------|
| May 1 | Design & implement 06b_immune_subclustering.py | 3 | Script + docstring |
| May 2 | Implement Snakemake rule + config integration | 2 | rules/06b_*.smk + config updates |
| May 3 | Test on lung dataset (268k cells, 15 clusters) | 2 | Validation report, purity metrics |
| May 4 | Refine markers based on lung results; update documentation | 2 | Enhanced marker panels, IMMUNOLOGY_GUIDE.md |
| May 5 | Finalize + commit to repo | 1 | Clean git history, ready for liver |

**Total effort:** ~10 hours (includes buffer for debugging)

---

## Part 4: Scalability for Future Immune Panels

### 4.1 Configuration-Driven Approach

The implementation uses **zero hard-coded markers**. All cell type definitions live in `config.yaml`, enabling:

```yaml
# For future immune panels: just update markers, no code changes
immune_markers:
  # CD8 T cells with exhaustion markers (5K panel has PDCD1, HAVCR2, etc.)
  CD8_T_cell_exhausted: ["CD8A", "PDCD1", "HAVCR2", "TIGIT", "LAG3"]
  CD8_T_cell_activated: ["CD8A", "IFNG", "GZMA", "GZMB", "IL2"]
  # Macrophage metabolic states (5K panel has LDHA, IDO1, etc.)
  M1_Macrophage: ["CD68", "LDHA", "PFKFB3", "NOS2", "TNF"]
  M2_Macrophage: ["CD68", "IDO1", "IL4RA", "ARG1", "IL10"]
```

### 4.2 Hepatology-Specific Adaptations

For liver immune panels (future), markers would include:

```yaml
# Liver Macrophage Subtypes
Kupffer_cell:
  markers: ["CD68", "MARCO", "CLEC4F", "VSIG4"]
  tissue: "liver"
  
Portal_macrophage:
  markers: ["CD68", "ITGAM", "CCL2", "MRC1"]
  tissue: "liver"
  
Stellate_cell:
  markers: ["COL1A1", "ACTA2", "DES", "PDGFRA"]
  tissue: "liver"
```

---

## Part 5: Documentation for Future Users

**New file:** `docs/IMMUNE_ANNOTATION_GUIDE.md` (to be written)
- How to interpret immune purity reports
- How to customize markers for new panels
- Troubleshooting guide (e.g., "All cells scored as Treg—what's wrong?")
- Best practices for immune cell definition

---

## Success Criteria

By end of Week 2:

✅ Immune subclustering step integrated into Snakemake DAG  
✅ Lung dataset: immune purity improved to ≥70% (from 15-59%)  
✅ Full validation report generated (heatmap + metrics)  
✅ Configuration-driven markers allow future customization  
✅ Documentation enables immune panel projects  
✅ All code committed to main branch  

---

## References & Literature

1. **T cell markers:**
   - Zheng et al. (2017). "Massively parallel digital transcriptional profiling of single cells." Nature Communications.
   - Human Protein Atlas immune markers.

2. **Macrophage polarization:**
   - Murray et al. (2014). "Macrophage Activation and Polarization." Nature Reviews Immunology.
   - scRNA-seq consensus on M1/M2 signatures.

3. **Xenium-specific considerations:**
   - 10x Genomics Xenium validation studies (preprints/pubs).
   - Spatial transcriptomics immune cell analysis (Weber et al., Nat Rev Methods).

---

**Document prepared by:** Claude Code  
**For:** Miguel Ángel Díaz-Campos (INMEGEN)  
**Next review:** After Week 2 implementation
