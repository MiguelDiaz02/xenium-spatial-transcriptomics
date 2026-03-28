# Xenium Spatial Transcriptomics Project

A curated knowledge base and analysis resource for **10x Genomics Xenium** spatial transcriptomics (ST) data, with a focus on liver and lung tissue applications.

## Overview

This repository consolidates:
- A reproducible Python analysis pipeline (`recode_st`) for Xenium 5k data
- Reference codebases from the spatial transcriptomics community (as git submodules)
- Relevant bibliography (PDFs)
- Custom gene panel proposals for lung and liver Xenium experiments
- Visual reference materials and presentation slides

## Repository Structure

```
Xenium_project/
├── xenium_code_knowledge/          # Reference codebases (git submodules)
│   ├── Xenium_5k_analysis_pipeline/    # Primary pipeline: recode_st package
│   ├── sopa/                           # Technology-agnostic spatial omics pipeline
│   ├── spatialdata_xenium_explorer/    # SpatialData → Xenium Explorer converter
│   ├── Xenium_benchmarking/            # Benchmarking 25+ Xenium datasets
│   ├── XeniumIO/                       # Bioconductor R package for Xenium I/O
│   ├── XeniumSpatialAnalysis/          # Spatial gradient & clustering (R)
│   ├── liver_ped_map/                  # Pediatric liver atlas + IFALD spatial analysis
│   ├── analysis_guides/                # Official 10x Genomics tutorial notebooks
│   └── spatialdata-notebooks/          # SpatialData community notebooks
│
├── xenium_bibliography_knowledge/  # Key scientific papers (PDF)
├── my_xenium_panel_markers/        # Custom Xenium gene panel proposals
│   ├── Lung panel (immunology focus)
│   └── Liver panel (August 2025)
├── my_project_objectives/          # Project proposal and objectives
├── slides_visual_references/       # Reference presentations and slides
└── how_to_create_slides/           # Scientific presentation design guides
```

## Primary Pipeline: `recode_st`

Located in `xenium_code_knowledge/Xenium_5k_analysis_pipeline/`, this is a fully reproducible Python package for end-to-end Xenium analysis.

### Key features
- **SpatialData** (Zarr-backed) as the core data format throughout
- Sequential pipeline orchestrated via `__main__.py`
- Typed configuration via Pydantic (`config.py`)
- Integration with **scVI-tools**, **Squidpy**, and **MuSpAn** for spatial statistics
- Cell-type annotation, doublet detection, denoising (ResolVI), and pseudobulk analysis

### Pipeline stages

| Stage | Module(s) | Description |
|-------|-----------|-------------|
| Ingest | `format_data.py`, `subsample_data.py` | Xenium output → SpatialData Zarr |
| QC | `qc.py` | Quality metrics and filtering |
| Dimensionality reduction | `dimension_reduction.py` | PCA / UMAP / Leiden clustering |
| Annotation | `annotate.py` | Cell-type labeling |
| Integration | `integrate_scvi.py`, `integrate_ingest.py` | scVI and scanpy ingest |
| Spatial analysis | `spatial_statistics.py`, `muspan.py` | Squidpy + MuSpAn |
| Visualization | `view_images.py` | Spatial image overlays |
| Downstream | `drug2cell.py`, `doublet_identification.py`, `psuedobulk.py`, `denoise_resolvi.py` | Downstream analyses |

### Setup

```bash
cd xenium_code_knowledge/Xenium_5k_analysis_pipeline

# Conda (recommended)
conda env create -f environment.yml
conda activate xenium_5k_venv
pip install --no-build-isolation --no-deps -e .

# Or venv
python -m venv recode_st && source recode_st/bin/activate
pip install -r requirements.txt && pip install -e .

# Dev dependencies
make install-dev
```

### Common commands

```bash
# Run tests
pytest tests/
pytest -v -p no:warnings --cov=src --cov-report=html --doctest-modules

# Lint
ruff check src/
ruff format src/

# Docs server (mkdocs at localhost:8000)
make serve
```

## Reference Codebases

All submodules in `xenium_code_knowledge/` are third-party repositories. Clone with submodules:

```bash
git clone --recurse-submodules https://github.com/MiguelDiaz02/<repo-name>.git
# or, if already cloned:
git submodule update --init --recursive
```

| Submodule | Language | Purpose |
|-----------|----------|---------|
| `Xenium_5k_analysis_pipeline` | Python | Primary analysis pipeline (`recode_st`) |
| `sopa` | Python | Technology-agnostic ST pipeline with Snakemake |
| `spatialdata_xenium_explorer` | Python | Convert SpatialData → Xenium Explorer format |
| `Xenium_benchmarking` | Python | Benchmarking across 25+ datasets; `xb` module |
| `XeniumIO` | R | Bioconductor package for reading Xenium data |
| `XeniumSpatialAnalysis` | R | Spatial gradient and clustering analysis |
| `liver_ped_map` | Python/R | Pediatric liver atlas and IFALD spatial analysis |
| `analysis_guides` | Python/R | Official 10x Genomics tutorial notebooks |
| `spatialdata-notebooks` | Python | Community notebooks for the SpatialData ecosystem |

## Gene Panel Proposals

Custom Xenium gene panel designs are in `my_xenium_panel_markers/`:

- **Lung panel** (`Propuesta panel Pulmón 2025`): Immunology-focused panel including immune cell markers, fibrotic lung markers, and airway/alveolar cell type markers.
- **Liver panel** (`Propuesta panel Xenium Liver August 2025`): Liver cell type markers for spatial transcriptomics profiling.

## Key Dependencies

- [`spatialdata`](https://github.com/scverse/spatialdata) / [`spatialdata-io`](https://github.com/scverse/spatialdata-io) — core spatial data format
- [`scanpy`](https://scanpy.readthedocs.io/) + [`anndata`](https://anndata.readthedocs.io/) — single-cell analysis
- [`squidpy`](https://squidpy.readthedocs.io/) — spatial statistics
- [`scvi-tools`](https://scvi-tools.org/) — probabilistic models for single-cell omics
- [`MuSpAn`](https://muspan.gitlab.io/) — multiscale spatial analysis
- `torch` — required by scVI and ResolVI

## Bibliography

Key references are stored as PDFs in `xenium_bibliography_knowledge/`. Topics covered include:
- Spatial transcriptomics benchmarking and methodology
- Pediatric liver atlas and spatial profiling
- Xenium platform validation and analysis workflows
- Cell-type deconvolution and integration methods

## License

This repository aggregates resources from multiple sources. Each submodule retains its original license. Please refer to individual submodule directories for their respective licenses.
