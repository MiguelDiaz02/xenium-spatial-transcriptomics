# Immune Cell Annotation Guide

**For:** Researchers using Xenium spatial transcriptomics data  
**Updated:** 2026-04-24  
**Version:** 1.0

---

## Quick Start

The Xenium pipeline includes a **step 06b: Immune Cell Subclustering** that automatically improves immune cell type annotation from generic labels (e.g., "T_cell") to granular subtypes (e.g., "CD8_T_cell", "Treg", "M1_Macrophage").

### To enable immune subclustering:

Ensure your `config.yaml` has:

```yaml
immune_annotation:
  enabled: true
  immune_cell_types: ["T_cell", "B_cell", "Macrophage", "NK_cell", "Mast_cell"]
  leiden_resolution_immune: 1.2
  validation_threshold: 0.7
```

### Run the pipeline:

```bash
cd pipeline
snakemake --configfile config/config_lung.yaml --cores 8
# Or for liver:
snakemake --configfile config/config_liver.yaml --cores 8
```

The step 06b runs automatically after annotation (step 06) and produces:
- Enhanced AnnData with `obs["cell_type_immune_granular"]` column
- Validation report: `results/06b_immune_subclustering_report.html`

---

## Understanding the Outputs

### Main Output Column: `cell_type_immune_granular`

This column contains fine-grained immune cell type labels:

```
CD8_T_cell        → Cytotoxic T cells
CD4_T_cell        → Helper T cells
Treg              → Regulatory T cells (FOXP3+)
M1_Macrophage     → Pro-inflammatory macrophage
M2_Macrophage     → Anti-inflammatory macrophage
Kupffer_cell      → Resident liver macrophage (liver panel)
NK_cell           → Natural killer cells
B_cell            → B lymphocytes
Plasma_cell       → Antibody-secreting cells
Dendritic_cell    → Antigen-presenting cells
Non-immune        → Non-immune cell types (epithelial, etc.)
```

### Supporting Columns

**`immune_purity`** (0-1 range)
- Quantitative score of how "pure" the annotation is
- Computed as: `max_score - median_other_scores`
- **0.9-1.0:** Excellent — cell clearly matches its assigned type
- **0.7-0.9:** Good — confident assignment
- **0.5-0.7:** Ambiguous — cell may express markers of multiple types (REVIEW)
- **<0.5:** Poor — uncertain assignment (FAIL)

**`immune_confidence`** (categorical)
- **PASS:** Purity ≥ 0.7 — high confidence annotation
- **REVIEW:** Purity 0.5-0.7 — ambiguous, warrants inspection
- **FAIL:** Purity < 0.5 — low confidence, may need remerging
- **N/A:** Non-immune cells (not in immune subset)

**`leiden_immune`** (categorical)
- Leiden cluster ID from immune-specific reclustering
- Useful for identifying co-clusters with similar transcriptomes

---

## Interpreting the Validation Report

### Report Contents

The HTML report (`06b_immune_subclustering_report.html`) includes:

**1. Immune distribution bar chart**
```
Shows count of each granular immune type:
  CD8_T_cell:   15,234 cells
  CD4_T_cell:    8,910 cells
  M1_Macrophage: 3,421 cells
  ...
```

**2. Purity histogram**
```
Distribution of immune_purity scores across all immune cells.
- Vertical green line: validation_threshold (default 0.7)
- Vertical orange line: REVIEW threshold (0.5)
- A right-skewed distribution (most cells >0.7) → good annotation
- A bimodal/flat distribution → marker panels may need refinement
```

**3. Purity by cluster boxplot**
```
Shows median purity per leiden_immune cluster.
- All boxes above green line → all subclusters are confident
- Some boxes below green line → those subclusters may need re-annotation
```

**4. Confidence summary bar chart**
```
Counts of PASS / REVIEW / FAIL cells.
- Ideal: >80% PASS, <10% FAIL
- If >30% FAIL → marker panels may not fit your dataset
```

---

## Troubleshooting Guide

### Problem 1: "All cells scored as X_cell (100% assigned to one type)"

**Symptom:** Every immune cell gets the same label (e.g., "Macrophage").

**Cause:** 
- Marker genes missing from your panel
- Marker genes don't discriminate between types
- Panel design issue (e.g., no T-cell markers)

**Solution:**
1. Check the script log: which marker genes were found vs. missing?
2. Update `immune_markers` in config to use genes present in your panel
3. Example: If your 289-gene panel lacks CD4/CD8/FOXP3, use alternative markers:
   ```yaml
   CD8_T_cell: [CD8A, GZMA, PRF1]      # Use effector markers only
   CD4_T_cell: [CD4, IL7R, SELL]        # Use alternative CD4 markers
   ```

### Problem 2: "All cells have purity < 0.5 (FAIL)"

**Symptom:** Validation histogram shows purity heavily left-shifted; >80% FAIL.

**Cause:**
- Immune populations are highly mixed / don't have distinct expression
- Marker panels are conflicting (cells express markers for multiple types)
- Dataset is unusual (e.g., immune-rich tissue; unusual cell activation state)

**Solution:**
1. Inspect the validation heatmap (in PDF): Which markers are being expressed broadly?
2. Visually check UMAP (immune subset): Are there clear clusters?
3. Option A: Relax validation threshold
   ```yaml
   validation_threshold: 0.5  # REVIEW replaces PASS at 0.5
   ```
4. Option B: Simplify marker panels (use fewer, more specific markers)
5. Option C: Don't use step 06b; rely on global `cell_type` instead

### Problem 3: "One specific immune type (e.g., Treg) has 0 cells"

**Symptom:** A valid immune subtype is never assigned.

**Cause:**
- Marker genes for that type are missing from panel
- That type doesn't exist in your tissue (expected)

**Solution:**
1. Check if Treg markers (FOXP3, IL2RA) are in your panel
2. If missing, remove that subtype from config:
   ```yaml
   immune_markers:
     CD8_T_cell: [...]
     # Treg removed (no FOXP3 in panel)
     B_cell: [...]
   ```
3. If present and still 0 cells → that cell type may not exist in your tissue (expected)

### Problem 4: "Validation report says REVIEW / FAIL — what do I do?"

**Symptom:** 30-50% of immune cells have `immune_confidence == 'REVIEW'` or `'FAIL'`.

**Cause:**
- Immune subtypes are genuinely mixed in your tissue
- OR marker panels are overlapping / not discriminative enough

**Solution — Strategy A: Accept ambiguity**
```python
# In downstream analysis, treat REVIEW/FAIL cells as "Unknown_immune"
adata.obs['cell_type_for_analysis'] = adata.obs['cell_type_immune_granular']
adata.obs.loc[adata.obs['immune_confidence'] != 'PASS', 'cell_type_for_analysis'] = 'Unknown_immune'
```

**Solution — Strategy B: Refine markers**
1. Manually inspect UMAP colored by `leiden_immune`
2. For a low-purity cluster, compute top differential genes
3. Update that subtype's marker panel with those genes
4. Re-run pipeline

---

## Customizing for Your Panel

### Step 1: Identify Available Markers

```python
import anndata as ad

adata = ad.read_h5ad('your_data.h5ad')
print(adata.var_names)  # List all genes in panel
```

### Step 2: Design Your Immune Subtypes

Example: If you have a **lung panel with immune focus**, you might define:

```yaml
immune_markers:
  # T cells (CD3, CD4, CD8 available in lung panel)
  CD8_T_cell_cytotoxic: [CD8A, CD8B, GZMA, GZMB, PRF1, GNLY]
  CD8_T_cell_exhausted: [CD8A, PDCD1, HAVCR2, TIGIT, LAG3]  # Add exhaustion markers
  CD4_T_cell_naive: [CD4, IL7R, CCR7, SELL]
  CD4_T_cell_activated: [CD4, IFNG, IL2, TNF]
  Treg: [FOXP3, IL2RA, CTLA4, IKZF2]
  
  # Macrophages (detailed polarization for lung immune panels)
  M1_Macrophage: [CD68, LDHA, PFKFB3, NOS2, TNF, IL6]
  M2_Macrophage: [CD68, IDO1, IL4RA, ARG1, IL10, TGM2]
  Alveolar_macrophage: [CD68, MARCO, MRC1, MERTK]  # Tissue-resident
  
  # Other
  NK_cell_activated: [GNLY, NKG7, IFNG, GZMA]
  B_cell: [MS4A1, CD19, IGHM]
```

### Step 3: Update Config & Re-run

```yaml
# config/config_mylung.yaml
immune_annotation:
  enabled: true
  immune_cell_types: ["T_cell", "B_cell", "Macrophage", "NK_cell"]
  immune_markers:
    CD8_T_cell_cytotoxic: [...]  # Your custom markers
    # ... rest
```

```bash
snakemake --configfile config/config_mylung.yaml --cores 8
```

---

## Adapting to Liver (or Other Tissues)

### Pre-made Configuration

We provide `config/config_liver.yaml` with liver-specific immune types:

```yaml
immune_markers:
  Kupffer_cell: [CD68, MARCO, CLEC4F, VSIG4, CD163]  # Resident liver Mφ
  Portal_macrophage: [CD68, ITGAM, CCL2, MRC1]       # Infiltrating Mφ
  CD8_T_cell: [CD8A, CD8B, GZMA, GZMB, PDCD1]
  Treg: [FOXP3, IL2RA, CTLA4, IKZF2]
  # ... etc
```

### For Other Tissues

Create a new config following this template:

```yaml
# config/config_pancreas.yaml (example)
immune_annotation:
  enabled: true
  immune_cell_types: ["T_cell", "B_cell", "Macrophage", "Dendritic_cell"]
  leiden_resolution_immune: 1.2
  validation_threshold: 0.7
  immune_markers:
    # Pancreatic-specific immune subtypes
    Pancreatic_stellate_cell_associated_macrophage: [CD68, TPSAB1, ...]
    Activated_T_cell: [CD3D, IFNG, GZMA, ...]
    # ... rest
```

---

## FAQ

**Q: Does step 06b require GPU?**  
A: No, it runs on CPU (threads-based parallelization). GPU is only used for upstream steps (segmentation, ResolVI).

**Q: What's the computational cost?**  
A: For 268k cells (lung dataset): ~3-5 minutes on 8 threads. Scales ~linearly with cell count.

**Q: Can I disable step 06b?**  
A: Yes:
```yaml
immune_annotation:
  enabled: false
```
The global `cell_type` from step 06 will be used for downstream analysis.

**Q: Can I modify immune markers without re-running the whole pipeline?**  
A: No, you must re-run step 06b. Use:
```bash
snakemake --configfile config/config.yaml --cores 8 --forcerun 06b_immune_subclustering
```

**Q: How do I combine global and granular annotations in my analysis?**  
A: Use both columns:
```python
# For spatial analysis that needs broad categories
adata.obs['cell_type_broad'] = adata.obs['cell_type']  # global

# For within-immune analysis that needs granularity
immune_only = adata[adata.obs['immune_confidence'] != 'N/A']
immune_only.obs['subtype'] = immune_only.obs['cell_type_immune_granular']
```

---

## Citation & Reproducibility

To cite this immune annotation workflow, reference:

**XENIUM_ANNOTATION_STRATEGY.md** (v1.0, 2026-04-24)  
Miguel Ángel Díaz-Campos, INP (Instituto Nacional de Pediatría)

Include the config file and immune markers used in your Methods:

> "Immune cells were subsetted and re-clustered at resolution 1.2, then assigned granular subtypes using marker scoring (Supplementary Table: immune_markers). Purity was validated using a 0.7 threshold; cells below 0.5 were flagged as ambiguous."

---

## Support & Feedback

For issues, questions, or marker suggestions:
- Check this guide and the troubleshooting section
- Review `IMMUNE_ANNOTATION_STRATEGY.md` for technical details
- Inspect validation heatmaps in `06b_immune_subclustering_report.pdf`

---

**Last Updated:** 2026-04-24  
**Maintainer:** Claude Code (xenium-pipeline)
