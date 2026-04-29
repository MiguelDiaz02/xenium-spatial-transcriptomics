# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a collection of spatial transcriptomics (ST) analysis resources focused on 10x Genomics Xenium data. There are two primary codebases:

1. **`xenium_code_knowledge/Xenium_5k_analysis_pipeline/`** — Reproducible Python pipeline (`recode_st` package) for modular, reusable ST analysis
2. **`proyecto_demo_xenium/`** — Snakemake-based pipeline for end-to-end analysis (human lung cancer pilot, scales to liver via config)

Supporting resources: reference codebases (git submodules), scientific papers, gene panel proposals, and presentation materials.

## Primary Project: recode_st (Python Package)

Located in `xenium_code_knowledge/Xenium_5k_analysis_pipeline/`. This is a modular, importable package for ST analysis with full test coverage and documentation.

### Setup

```bash
cd xenium_code_knowledge/Xenium_5k_analysis_pipeline

# Conda (recommended; Python 3.10+)
conda env create -f environment.yml
conda activate xenium_5k_venv
pip install --no-build-isolation --no-deps -e .

# Or venv
python -m venv recode_st && source recode_st/bin/activate
pip install -r requirements.txt && pip install -e .

# Dev dependencies + docs
pip install -e .[dev]
# or: make install-dev
```

### Development Workflow

```bash
# Pre-commit hooks (runs ruff, format, markdown lint, YAML format on commit)
pre-commit install
# To bypass in a pinch: git commit --no-verify (not recommended)

# Tests (must activate venv first)
pytest tests/
pytest tests/test_config.py                        # Single file
pytest -v -p no:warnings --cov=src --cov-report=html --doctest-modules

# Lint & format (ruff; configured in pyproject.toml)
ruff check src/
ruff format src/
ruff check --fix src/      # Auto-fix imports, upgrades, etc.

# Docs server (mkdocs at localhost:8000)
make serve
```

### Architecture

Pipeline stages in `src/recode_st/`:

1. **Entry & Config** — `__main__.py` orchestrates; `config.py` uses Pydantic for typed configuration
2. **Ingest** — `format_data.py` (Xenium output → SpatialData Zarr); `subsample_data.py` for downsampling
3. **QC** — `qc.py` computes and filters quality metrics
4. **Dimensionality reduction & clustering** — `dimension_reduction.py` (PCA/UMAP/Leiden)
5. **Annotation** — `annotate.py` for cell-type labeling
6. **Integration** — `integrate_scvi.py` (scVI-tools), `integrate_ingest.py` (scanpy ingest)
7. **Spatial analysis** — `spatial_statistics.py` (Squidpy), `muspan.py` / `ms_spatial_graph.py` / `ms_spatial_stat.py` (MuSpAn, optional)
8. **Visualization** — `view_images.py` for spatial image overlays
9. **Downstream** — `drug2cell.py`, `doublet_identification.py`, `pseudobulk.py`, `denoise_resolvi.py`

**Data format:** All intermediate and final results are **SpatialData** objects (Zarr-backed).

### Key Dependencies

- `spatialdata` / `spatialdata-io` — core spatial data format
- `scanpy` + `anndata` — single-cell analysis
- `squidpy` — spatial statistics
- `torch` — required by scVI-tools and ResolVI
- `MuSpAn` — optional; for multiscale spatial analysis (see below)

### MuSpAn (Optional Spatial Graph Module)

To enable MuSpAn-based spatial analysis:

```bash
pip install -e ".[muspan]"
# or manually:
pip install -r requirements-muspan.txt
```

MuSpAn installation is non-standard; if it fails, debug with:
```bash
pip install --upgrade pip
pip install --no-cache-dir -e ".[muspan]"
```

## Secondary Project: Snakemake Pipeline

Located in `proyecto_demo_xenium/`. A modular, reproducible Snakemake workflow for end-to-end Xenium analysis. Designed as a pilot on human lung cancer (289-gene panel) but scales to liver (5K panel) via config only.

See `proyecto_demo_xenium/CLAUDE.md` for detailed pipeline documentation, including setup, running, adapting to liver data, and layer conventions.

Quick reference:

```bash
cd proyecto_demo_xenium/pipeline
conda env create -f envs/xenium_pipeline.yaml
conda activate xenium_pipeline

# Dry run — check DAG without executing
snakemake --configfile config/config_lung.yaml -n

# Full run (8 cores)
snakemake --configfile config/config_lung.yaml --cores 8

# Resume from a specific step
snakemake --configfile config/config_lung.yaml --cores 8 --forcerun qc
```

### Pipeline Status (2026-04-23)

**All 7 essential steps completed successfully.** Human lung cancer (V1 panel) pilot analyzed: 268,034 cells, 289 genes, 15 Leiden clusters, 7 annotated cell types.

**Critical fixes applied & validated:**

1. **spatialdata-io ≥0.1.6** — Fixed `imread() TypeError` from obsolete imageio API
2. **Segmentation: Cellpose → xoa** — Resolved ViT-SAM incompatibility with DAPI uint16 imagery
3. **SpatialData API v2 compliance** — Updated all 7 downstream scripts (sdata.table → sdata.tables[dict])
4. **Annotation: CellTypist → manual markers** — Model unavailability resolved via predefined gene signatures
5. **Categorical dtype enforcement** — Squidpy compatibility for cell_type columns

**Output:** `human_lung_cancer/results/sdata.zarr` (6.2 GB) — complete SpatialData object with UMAP, Leiden clusters, 7 cell type annotations, ResolVI-denoised expression, and raw counts preserved.

**Steps 08-09 status:** Disabled (optional, complex dependencies). Can be enabled post-debugging if needed.

For detailed completion report, see `XENIUM_PIPELINE_COMPLETION_REPORT.md`.

## Presentation Scripts (Root Level)

Python scripts for generating Xenium presentation materials. All require a virtual environment with `python-pptx`, `python-docx`, `cairosvg`, `requests`, and `lxml`.

```bash
# Set up a simple venv for presentations (if needed)
python -m venv pptx_env && source pptx_env/bin/activate
pip install python-pptx python-docx cairosvg requests lxml

# Generate slides (version 1 — dark navy theme, 3 slides)
python create_xenium_slides.py

# Generate slides (version 2 — redesigned, FlashTalk-MADC inspired)
python create_xenium_slides_v2.py

# Generate speaker notes / script
python create_script_docx.py

# Merge multiple PPTX files
python merge_pptx.py
```

| Script | Output |
|--------|--------|
| `create_xenium_slides.py` | `xenium_presentacion_final.pptx`, `xenium_presentacion_grupo_*.pptx` |
| `create_xenium_slides_v2.py` | Redesigned PPTX (columns, pill labels, banners) |
| `create_script_docx.py` | `guion_xenium.docx` (speaker notes) |
| `merge_pptx.py` | Combined slides (requires multiple PPTX sources) |

Slide preview PNGs are cached in `slide_previews/`.

## Git Submodules (Reference Codebases)

All under `xenium_code_knowledge/` as git submodules. Initialize and update with:

```bash
# Clone with submodules (first time)
git clone --recurse-submodules <repo-url>

# Or update existing repo
git submodule update --init --recursive

# Check status
git config --file .gitmodules --name-only --get-regexp path
```

| Directory | Language | Purpose |
|-----------|----------|---------|
| `sopa/` | Python | Technology-agnostic spatial omics pipeline with CLI and Snakemake |
| `spatialdata_xenium_explorer/` | Python | Convert SpatialData → Xenium Explorer format |
| `Xenium_benchmarking/` | Python | Benchmarking across 25+ datasets; `xb` module |
| `XeniumIO/` | R | Bioconductor package for reading Xenium data |
| `XeniumSpatialAnalysis/` | R | Spatial gradient and clustering analysis |
| `liver_ped_map/` | Python/R | Pediatric liver atlas + IFALD spatial analysis |
| `analysis_guides/` | Python/R | Official 10x Genomics tutorial notebooks |
| `spatialdata-notebooks/` | Python | SpatialData community notebooks |

## Non-code Resources

| Path | Content |
|------|---------|
| `xenium_bibliography_knowledge/` | Key papers as PDFs (ST benchmarking, liver atlas, Xenium validation, deconvolution) |
| `my_xenium_panel_markers/` | Custom Xenium gene panel proposals: lung (immunology-focused) and liver (August 2025) |
| `my_project_objectives/` | Project proposal PDF (`Anteproyecto_MADC-1.pdf`) |
| `slides_visual_references/` | Reference PPTX and PDF presentations |
| `how_to_create_slides/` | Presentation design guides |

## Week 3 Status: Phase 1A & 1B COMPLETE ✅

**Date**: 2026-04-24  
**Completion**: Spatial-Aware DGE Analysis + Within-Celltype Validation

### Phase 1A Results
- **288/289 genes (99.7%)** with significant spatial autocorrelation (Moran's I, p<0.05)
- **Wilcoxon DE analysis** for all 7 cell types (B_cell, Endothelial, Epithelial, Macrophage, NK_cell, T_cell, Tumor)
- **100% spatial enrichment** in top-50 DE markers per cell type
- **Outputs**: 12 CSV files + 4 publication figures (PNG) + 1 comprehensive PDF with captions
- **Location**: `human_lung_cancer/results/02_biology/immune_DE/` + `figures/phase1a_spatial_de/`

### Phase 1B Results
- **325/350 robust markers (92.86%)** validated by DE + within-celltype spatial signal
- **Phase 1A ↔ 1B concordance**: r>0.85 (highly consistent, biologically valid)
- **54 genes flagged** as scattered expression (non-clustered, deprioritized)
- **Outputs**: 4 CSV files + 3 publication figures (PNG) + 1 comprehensive PDF with captions
- **Location**: `human_lung_cancer/results/02_biology/immune_DE/phase1b_validation/` + `figures/phase1b_spatial_de_validation/`

### Next Phase: Phase 2 (CCC Analysis)
- **Method**: DeepTalk (graph attention networks)
- **Focus**: Ligand-receptor interactions, M2 macrophage hub validation
- **Timeline**: 4-6 hours
- **Status**: Ready to begin

### Key Scripts Created
- `pipeline/scripts/analysis/week3_01a_dge_squidpy.py` — Phase 1A spatial-aware DE
- `pipeline/scripts/analysis/week3_01a_visualizations.py` — Phase 1A figures
- `pipeline/scripts/analysis/week3_01b_spatial_de_validation.py` — Phase 1B validation
- `pipeline/scripts/analysis/week3_01b_visualizations.py` — Phase 1B figures

### File Organization
```
proyecto_demo_xenium/
├── human_lung_cancer/results/02_biology/immune_DE/          [16 CSVs]
├── human_lung_cancer/results/figures/
│   ├── phase1a_spatial_de/                                   [4 PNGs + PDF]
│   ├── phase1b_spatial_de_validation/                        [3 PNGs + PDF]
│   ├── phase2_ccc/                                           [empty - next]
│   ├── phase3_spatial/                                       [empty - next]
│   └── phase4_summary/                                       [empty - next]
├── pipeline/scripts/analysis/week3_*.py                      [4 scripts]
├── WEEK3_FILE_LOCATIONS.txt                                  [quick reference]
└── WEEK3_ANALYSIS_SUMMARY.md                                 [detailed summary]
```

## Autonomous Operation Guidelines

When operating autonomously (sol rojo mode):
- **Data format:** SpatialData (Zarr). Never convert to other formats unless explicitly requested.
- **Environment:** Activate `xenium_5k_venv` (recode_st) or `xenium_pipeline` (Snakemake) before running code.
- **Data safety:** No destructive ops on raw data — write all outputs to `results/` or a dedicated output directory.
- **Failures:** Try one alternative; log errors to `logs/` and continue.
- **Patterns:** Prefer patterns from `src/recode_st/` (Python) or `pipeline/scripts/` (Snakemake).
- **Completion:** Write summary to `logs/session_summary_YYYYMMDD.txt`.
- **Figure organization:** Always organize figures in `results/figures/` with subfolders per analysis phase
- **Documentation:** Update CLAUDE.md and memory files at end of session to preserve context
- **Reproducibility:** All scripts should be self-contained and include logging via `utils.logging.get_logger()`
