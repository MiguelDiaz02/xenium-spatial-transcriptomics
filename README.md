# Xenium Spatial Transcriptomics: Multi-Tool Consensus Pipeline

**A modular, reproducible Snakemake pipeline for 10x Genomics Xenium spatial transcriptomics analysis.** Designed for publication in *Nature Communications* with emphasis on triangulation and multi-tool consensus across 8 analytical families.

**Dataset:** 268,034 cells × 289 genes | Human lung cancer, Xenium v1 (pilot)  
**Target:** Scalable to Xenium Prime 5K liver data (config-only changes)  
**Hardware:** GPU-enabled (RTX 4500 Ada, 18GB VRAM recommended; CPU fallback available)

---

## Overview

This repository implements a **consensus-first philosophy** for spatial transcriptomics analysis, inspired by LIANA+ (aggregating multiple methods before ranking). Every biological finding is triangulated across ≥2 independent tools before reporting.

### 8-Family Analytical Framework

| Family | Focus | Method | Status |
|--------|-------|--------|--------|
| **F1** | Spatial domain detection | Banksy clustering + Leiden | ✅ Complete |
| **F2** | Spatially variable genes | Hotspot + nnSVG + Moran's I | 🔄 Scaffolded |
| **F3** | Cell-cell communication validation | LIANA+ bivariate + Spacia | ⏳ Pending |
| **F4** | Niche-dependent differential expression | NicheDE (Mason 2024) | ⏳ Pending |
| **F5** | Spatial pseudotime trajectories | Slingshot + tradeSeq + GAM | ⏳ Pending |
| **F6** | Cell type annotation x-validation | SingleR (optional) | ⏳ Pending |
| **F7** | Co-expression networks | ~~hdWGCNA~~ | ❌ Skipped (insufficient genes) |
| **F8** | Foundation model embeddings | Novae (Nat Methods 2025) | ✅ Complete |

---

## Project Status (2026-04-28)

### Completed ✅

- **Phase A**: Snakemake skeleton (steps 01–12) + hierarchical re-labeling (24 L2 cell types)
- **Phase B**: Spatial domains via F8 (Novae) + F1 (Banksy) with consensus (94.8k robust cells)
- **Phase 1A** (supplementary): DGEA benchmarked (288/289 genes spatially autocorrelated)
- **Phase 2B** (supplementary): CCC via LIANA+ (141 significant L-R interactions)

### Scaffolded 🔄

- **Phase C**: SVG consensus (`F2_hotspot_svg.py` ✅, `F2_nnsvg.R` ✅, `F2_svg_consensus.py` pending)

### Pending ⏳

- **Phase D**: Niche-DE (F4)
- **Phase E**: Spatial pseudotime (F5)
- **Phase F**: CCC validation (F3)
- **Phase G**: Annotation x-validation + master figure (F6)

---

## Installation

### Prerequisites

- Python ≥3.10, R ≥4.2
- Conda (Mambaforge recommended)
- CUDA ≥11.8 (GPU strongly recommended for Cellpose, ResolVI)
- 32+ GB RAM (64+ for reproducibility)

### Quick Start

```bash
# Clone the repository
git clone https://github.com/MiguelDiaz02/xenium-spatial-transcriptomics.git
cd xenium-spatial-transcriptomics/proyecto_demo_xenium

# Create and activate the conda environment
conda env create -f pipeline/envs/xenium_pipeline.yaml
conda activate xenium_pipeline

# For Phase C (SVG consensus), also create the R environment
conda env create -f pipeline/envs/xenium_R_analysis.yaml

# Verify GPU access
python -c "import torch; print('GPU available:', torch.cuda.is_available())"
```

#### Optional: Julia + Baysor (for transcript-level segmentation)

```bash
curl -sSL https://install.julialang.org | sh
julia -e 'using Pkg; Pkg.add("Baysor")'
pip install baysorpy
```

---

## Running the Pipeline

### Preprocessing (Steps 01–12)

```bash
cd proyecto_demo_xenium/pipeline

# Dry run (check DAG without execution)
snakemake --configfile config/config_lung.yaml -n

# Full run (8 cores, GPU auto-detected)
snakemake --configfile config/config_lung.yaml --cores 8

# Resume from a specific step (e.g., after fixing QC thresholds)
snakemake --configfile config/config_lung.yaml --cores 8 --forcerun qc

# Visualize DAG
snakemake --configfile config/config_lung.yaml --dag | dot -Tpdf > dag.pdf
```

### Analysis Phases (F1–F8)

Each phase has its own Snakemake rule file under `pipeline/analysis_rules/`:

```bash
# Phase B (Spatial Domains F1+F8) — already complete
snakemake --configfile config/config_lung.yaml --cores 8 -s pipeline/analysis_rules/F1_F8_spatial_domains.smk

# Phase C (SVG Consensus F2) — scaffolded, ready to execute
snakemake --configfile config/config_lung.yaml --cores 8 -s pipeline/analysis_rules/F2_svg_consensus.smk
```

---

## Data Format

### Input: Xenium Raw Output

- Directory structure: `experiment.xenium/` (ZIP or folder)
- Contains: morphology.ome.tif, cell_feature_matrix.h5, cell_boundaries.parquet, transcripts.parquet

### Processed Data: SpatialData Zarr

All intermediate and final results are stored in a **single SpatialData Zarr store**:

```
results/sdata.zarr/
├── tables/
│   ├── table_0.h5ad          # Main AnnData (268,034 cells × 289 genes)
│   │   ├── obs/              # Cell metadata
│   │   ├── X                 # Normalized + log1p expression
│   │   └── layers/
│   │       ├── "counts"      # Raw integer counts (for DE, scVI, ResolVI)
│   │       └── "denoised"    # ResolVI-denoised (visualization only)
├── shapes/                   # Cell boundaries (xoa segmentation)
├── images/                   # H&E, DAPI, morphology channels
└── labels/                   # Spatial domain masks (Phase B)
```

---

## Key Results

### Phase B: Spatial Domain Detection (✅ 2026-04-28)

- **F8 Novae:** 51.7s, 64d embeddings, hierarchical domains
- **F1 Banksy:** 531.9s, 14 spatial domains via Leiden
- **Consensus:** 94,799 robust cells (35.4% agreement), niche + transition zones
- **Output location:** `human_lung_cancer/results/03_phase3_spatial/`

### Phase 1A: Differential Gene Expression (Supplementary, ✅)

- **Result:** 288/289 genes (99.7%) with significant spatial autocorrelation (Moran's I)
- **Top markers per cell type** ranked by spatial signal
- **Output location:** `human_lung_cancer/results/02_biology/immune_DE_benchmarked/`

### Phase 2B: Cell-Cell Communication (Supplementary, ✅)

- **LIANA+ rank aggregation:** 141 significant L-R interactions
- **Key finding:** Weak checkpoint axis; M2 macrophages as suppression hub
- **Therapeutic implication:** M2 targeting + monocyte differentiation blocking likely more effective than anti-PD-1 monoterapy
- **Output location:** `human_lung_cancer/results/02_biology/ccc_liana/`

---

## Repository Structure

```
xenium-spatial-transcriptomics/
├── README.md                           # This file
├── proyecto_demo_xenium/               # Main Snakemake pipeline
│   ├── pipeline/
│   │   ├── Snakefile                   # Main orchestrator (steps 01–12)
│   │   ├── config/
│   │   │   ├── config_lung.yaml        # Active: 289-gene lung data
│   │   │   └── config_liver.yaml       # Future: 5K-gene liver data
│   │   ├── rules/                      # Preprocessing rules (01–09)
│   │   ├── analysis_rules/             # Analytical families (F1–F8)
│   │   ├── scripts/
│   │   │   ├── *_[01-09].py            # Preprocessing scripts
│   │   │   └── analysis/               # Phase 1A, 2B, B, C scripts
│   │   ├── utils/
│   │   │   ├── io.py                   # SpatialData helpers
│   │   │   └── logging.py              # Rich logging
│   │   └── envs/
│   │       ├── xenium_pipeline.yaml    # Main environment
│   │       ├── xenium_R_analysis.yaml  # R for Phases C–E
│   │       └── spacia.yaml             # Spacia (isolated)
│   ├── human_lung_cancer/
│   │   └── results/                    # All outputs (sdata.zarr, CSVs, figures)
│   ├── SENDA_DORADA.md                 # Golden path: project state & phases
│   └── CLAUDE.md                       # Developer notes & reproducibility
├── xenium_code_knowledge/              # Reference codebases (git submodules)
├── xenium_bibliography_knowledge/      # Key papers (PDFs)
└── my_xenium_panel_markers/            # Custom gene panel proposals
```

---

## Configuration & Adaptation

### For Lung Data (Current)

```yaml
# config/config_lung.yaml
data:
  xenium_dir: "human_lung_cancer/rawdata/"
  output_dir: "human_lung_cancer/results/"

annotation:
  celltypist_model: "Immune_All_Low"

segmentation:
  method: "auto"  # Auto-selects: Cellpose (current)
```

### For Liver Data (Future, Config-Only)

```yaml
# config/config_liver.yaml
data:
  xenium_dir: "liver_healthy/rawdata/"
  output_dir: "liver_healthy/results/"

annotation:
  celltypist_model: "Liver_Human_PIP"

segmentation:
  method: "auto"  # Auto-selects: Baysor (5K-gene density)
```

**No rule or script changes required.**

---

## Performance Notes

### Large Dataset Optimization (>100k cells)

- **Step 09 (Doublet detection):** PCA components reduced from 30→15 for 268k cells (saves ~30 min)
- **Timeout:** 3600s | Memory: 32 GB
- **Tuning for >200k cells:** Further reduce PCA (n_components=10)

### Step Execution Times (Lung Dataset, GPU)

| Step | Time |
|------|------|
| 01 Ingest | 45s |
| 02 Segmentation (Cellpose) | 8–12 min |
| 03–05 QC, Preprocess, Reduction | 20–25 min |
| 06–07 Annotation, Denoising | 50–70 min |
| 08–09 Spatial, Downstream | 50–100 min |
| **Total** | **~3–3.5 hours** |

Phase B (F1+F8): 10 min | Phase 1A/2B: <1 min each

---

## Reproducibility

### Environment Locking

All conda environments are locked to exact versions:

```bash
conda env create -f pipeline/envs/xenium_pipeline.yaml
conda activate xenium_pipeline
```

### Snakemake DAG Visualization

```bash
snakemake --configfile config/config_lung.yaml --dag | dot -Tpdf > dag.pdf
```

---

## Citation

```bibtex
@software{Diaz-Campos2026,
  author       = {Díaz-Campos, Miguel Ángel},
  title        = {Xenium Spatial Transcriptomics: Multi-Tool Consensus Pipeline},
  year         = {2026},
  url          = {https://github.com/MiguelDiaz02/xenium-spatial-transcriptomics},
  organization = {Instituto Nacional de Pediatría (INP)},
  keywords     = {spatial transcriptomics, xenium, consensus methods, snakemake}
}
```

Also cite key publications:
- **LIANA+**: Dimitrov et al. Nature Methods 21 (2024)
- **Novae**: LIN et al. Nature Methods 21 (2025)
- **Banksy**: Foroutan et al. Nat Methods 20 (2023)

---

## Contact

**Author:** Miguel Ángel Díaz-Campos (computational biologist)  
**Affiliation:** Instituto Nacional de Pediatría (INP), México  
**Email:** miguelangeld00@gmail.com

---

## License

MIT License. See individual submodule directories for their respective licenses.
