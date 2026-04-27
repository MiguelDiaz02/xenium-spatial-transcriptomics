# Design Spec: Phase 3 — Spatial Mapping & Neighborhood Enrichment

**Date:** 2026-04-27  
**Project:** Xenium Lung Cancer — Spatial Transcriptomics Pipeline  
**Dataset:** 268,034 cells × 289 genes × 7 cell types  
**Status:** Approved — Ready for implementation

---

## 1. Context & Narrative Position

Phase 3 completes the spatial argumentative arc of the analysis:

| Phase | Question | Method |
|-------|----------|--------|
| 1A | Which genes have spatial signal? | DGEA + Moran's I |
| 2B | Which cell interactions are spatially contextualized? | CCC Hybrid (spatial-weighted) |
| **3** | **Where in the tissue do these interactions occur?** | **Spatial mapping + enrichment** |

Phase 3 answers: "In what functional tissue context do our validated markers and interactions operate?"

**Key inputs from prior phases:**
- `results/02_biology/immune_DE_benchmarked/` — Phase 1A spatial marker rankings
- `results/02_biology/ccc_hybrid_method/ccc_pairwise_interactions.csv` — Phase 2B CCC scores

---

## 2. Architecture & Snakemake Integration

**New files:**
```
pipeline/
├── rules/
│   ├── phase3_spatial_mapping.smk    # Tasks 3.1–3.3 rules
│   └── phase3_visualization.smk      # Integration visualization rule
├── scripts/analysis/
│   ├── phase3_01_spatial_gradients.py
│   ├── phase3_02_enrichment_weighted.py
│   ├── phase3_03_tissue_regions.py
│   ├── phase3_napari_review.py        # Interactive validation (manual step)
│   └── phase3_04_integration_viz.py
└── config/
    └── config_lung.yaml               # Add phase3 section
```

**Output directory:** `human_lung_cancer/results/03_phase3_spatial/`

**Snakemake execution:**
```bash
snakemake --configfile config/config_lung.yaml --cores 8 -R phase3_gradients phase3_enrichment phase3_regions phase3_viz
```

**DAG:**
```
phase3_gradients   →  immune_gradients.csv           ┐
phase3_enrichment  →  enrichment_spatial_weighted.csv ┤─► phase3_viz → integrated_tissue_map.png
phase3_regions     →  region_assignments.csv          │
                           │                           │
                    [MANUAL STEP: Napari review]       │
                           │                           │
                   region_assignments_validated.csv ──►┘
phase3_markers     →  region_markers.csv              ┘
```

**Napari review is a manual checkpoint** — not an automatic Snakemake rule. Pipeline execution:
1. Run `phase3_regions` → produces `region_assignments.csv`
2. User manually runs `phase3_napari_review.py` → reviews H&E overlay → saves `region_assignments_validated.csv`
3. Resume Snakemake → `phase3_viz` detects validated file and proceeds

---

## 3. Task 3.1 — Spatial Gradients (Methods A + B)

**Objective:** Map immune infiltration from tumor core → edge → healthy stroma.

### Method A: Squidpy Spatial Markers + Euclidean Distance
- Identify tumor boundary using `Epithelial` cell density + `sq.gr.spatial_markers()`
- Calculate euclidean distance of each cell to nearest tumor boundary point
- Bin distances: 0–100μm, 100–200μm, 200–300μm, 300–400μm, 400–500μm, >500μm
- Compute immune cell density per bin per type (T_cell, Macrophage, NK_cell, B_cell)

### Method B: KDE + Smoothed Gradients
- Apply `scipy.stats.gaussian_kde` to immune cell coordinates
- Compute spatial gradient (∇ density) using `numpy.gradient`
- Output gradient magnitude + direction per cell position

### Output: `immune_gradients.csv`
```
cell_id, x, y, cell_type, distance_to_tumor_um, distance_bin,
T_cell_density, Macrophage_density, NK_density, B_density,
gradient_magnitude, gradient_direction_deg
```

---

## 4. Task 3.2 — Neighborhood Enrichment Spatial-Weighted (Option Y)

**Objective:** Identify co-localizing cell type pairs, weighted by spatial proximity.

### Algorithm
```
For each cell pair (i, j) within knn_k neighbors:
  d = euclidean_distance(i, j)
  weight = exp(-d / bandwidth_um)
  enrichment_score[type_i, type_j] += weight

Normalize against null distribution:
  Permute cell_type labels 1000x
  z_score = (observed - mean_null) / std_null
```

### Cross-validation with Phase 2B
- Where `enrichment_zscore(A,B) > 2` AND `ccc_score(A→B) > threshold` → `validated_interaction = True`
- This produces **cross-validated interactions**: spatially co-localizing AND functionally interacting

### Config parameters (`config_lung.yaml`)
```yaml
phase3:
  enrichment:
    bandwidth_um: 50
    n_permutations: 1000
    knn_k: 15
    min_cells_per_type: 50
    ccc_score_threshold: null  # computed at runtime as 75th percentile of Phase 2B ccc_score distribution
```

### Spatial graph construction
- `sq.gr.spatial_neighbors(adata, n_neighs=15, radius=200)`
- Custom gaussian weights applied over this graph

### Output: `enrichment_spatial_weighted.csv`
```
cell_type_A, cell_type_B, enrichment_zscore, bandwidth_um,
ccc_score_phase2b, validated_interaction
```

---

## 5. Task 3.3 — Tissue Region Classification (Methods 1 + 3)

**Objective:** Define functional tissue zones and assign each cell to a region.

### Method 1: Spatial K-means + Marker Enrichment
```python
X_combined = np.hstack([
    coords_normalized * 0.7,          # spatial position (weight 0.7)
    phase1a_spatial_scores * 0.3      # top Moran's I markers (weight 0.3)
])
kmeans = KMeans(n_clusters=5, random_state=42)
```

- `n_clusters` configurable in `config_lung.yaml` (default: 5)
- After clustering: DE analysis per region to define marker enrichment
- Connects Phase 1A: verify spatially-variable markers enrich in coherent regions

### Method 3: Manual Validation with H&E Overlay (Napari)
- Load H&E image from `sdata.zarr` (pre-registered)
- Overlay K-means region colors on Napari viewer
- User validates / adjusts region labels and names
- Output: `region_assignments_validated.csv` (manually reviewed version)

```python
# phase3_napari_review.py
import napari
viewer = napari.Viewer()
viewer.add_image(he_image, name='H&E')
viewer.add_points(coords, properties={'region': region_labels},
                  face_color='region', size=3)
napari.run()
```

### Config parameters
```yaml
phase3:
  regions:
    n_clusters: 5
    spatial_weight: 0.7
    expression_weight: 0.3
    top_markers_per_region: 50
```

### Outputs
- `region_assignments.csv` — `cell_id, x, y, region_kmeans, region_validated, cell_type`
- `region_markers.csv` — `region, gene, log2fc, moran_i, pval_adj` (top 50 per region)

---

## 6. Integration Visualization

**File:** `phase3_integrated_tissue_map.png` (300 DPI, 20×16 inches)

**Layout 2×2:**
```
┌─────────────────────┬─────────────────────┐
│  A. Tissue Regions  │  B. Immune Gradient  │
│  (K-means validated │  (KDE density curves │
│   overlay on H&E)   │   per immune type)   │
├─────────────────────┼─────────────────────┤
│  C. Enrichment      │  D. Cross-Validated  │
│  Heatmap (weighted  │  Interactions        │
│  Z-score, 7×7)      │  (Phase2B × Phase3)  │
└─────────────────────┴─────────────────────┘
```

- **Panel A:** Spatial tissue map, Crameri colormap per region
- **Panel B:** KDE density curves, 4 immune types vs distance-to-tumor
- **Panel C:** Enrichment Z-score heatmap (Crameri diverging, center=0)
- **Panel D:** Bubble plot — X=CCC_score, Y=enrichment_zscore, size=n_cells, color=dominant region

Panel D is the narrative integration: interactions validated by both CCC (Phase 2B) and spatial co-localization (Phase 3).

---

## 7. Complete Data Flow

```
sdata.zarr
    │
    ├──► phase3_01_spatial_gradients.py
    │         ├── Phase 1A immune markers
    │         └──► immune_gradients.csv
    │
    ├──► phase3_02_enrichment_weighted.py
    │         ├── Phase 2B CCC scores
    │         └──► enrichment_spatial_weighted.csv
    │
    ├──► phase3_03_tissue_regions.py
    │         ├── Phase 1A spatial scores
    │         └──► region_assignments.csv
    │                   │
    │             [Napari review - manual]
    │                   │
    │             region_assignments_validated.csv
    │             region_markers.csv
    │
    └──► phase3_04_integration_viz.py
              ├── All 5 CSVs above
              └──► phase3_integrated_tissue_map.png
```

---

## 8. Deliverables Summary

| File | Task | Size (est.) | Purpose |
|------|------|-------------|---------|
| `immune_gradients.csv` | 3.1 | ~5 MB | Gradient curves per immune type |
| `enrichment_spatial_weighted.csv` | 3.2 | ~50 KB | Weighted enrichment matrix |
| `region_assignments.csv` | 3.3 | ~10 MB | Per-cell region labels |
| `region_assignments_validated.csv` | 3.3 | ~10 MB | Napari-reviewed labels |
| `region_markers.csv` | 3.3 | ~100 KB | Top 50 markers per region |
| `phase3_integrated_tissue_map.png` | All | ~3 MB | Manuscript figure |

---

## 9. Conda Environment

All scripts run under `xenium_pipeline` conda env (already configured).  
Additional packages needed:
- `scipy` — already present (KDE, distance calculations)
- `sklearn` — already present (KMeans)
- `napari` — installed on system (Napari review script)

No new environment changes required.
