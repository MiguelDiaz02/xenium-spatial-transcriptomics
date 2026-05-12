"""Synthetic mini-TMA cohort for P1 smoke testing.

Generates two SpatialData-backed TMAs (LUNG_AF, LIVER_AF) populated with
realistic gene-symbol expression so that LIANA+ can find ligand-receptor
pairs and decoupler/PyDESeq2 can do DE without crashing on edge cases.

Design (smoke variant — large enough to exercise stats, small enough to
run fast):

    LUNG_AF TMA:
        4 donors × 2 cores = 8 cores × 250 cells = 2000 cells
        Donors: 2 TBD (fibrotic_TBD) + 1 NS + 1 NN (both 'control')
        Cell types per core: AT1, Fibroblast, Macrophage_M1

    LIVER_AF TMA:
        6 donors × 2 cores = 12 cores × 250 cells = 3000 cells
        Donors: 2 TBD + 2 Control + 2 AlcCirh (all 3 conditions covered)
        Cell types per core: Hepatocyte, HSC_activated, Endothelial

Gene panel (~50 genes; intentionally includes LIANA consensus L-R pairs
for the fibrosis/immune axes the TBDs cohort cares about):

    Fibrosis:    TGFB1 TGFBR1 TGFBR2 COL1A1 COL3A1 ACTA2 PDGFB PDGFRB FAP
    Immune:      CD80 CD86 CTLA4 CD68 CD163 IL6 IL6R CXCL12 CXCR4
    Epithelial:  CDH1 EGFR AGER PDPN MUC1 ADAM17
    Endothelial: VWF PLVAP LYVE1
    Liver-specific: ALB APOA1 LRAT HNF4A
    Pan:         B2M GAPDH ACTB MALAT1 + 13 fillers

The expression model:
    - Baseline Poisson(λ=2)
    - Per-cell-type markers Poisson(λ=12)
    - Per-condition perturbation on a handful of fibrosis genes
      (TGFB1, COL1A1, ACTA2 upregulated in fibrotic_TBD)

This file is *only* for tests/smoke_test_p1.py. Production cohort YAML
is untouched (pipeline/config/cohort_TBDs.yaml).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.sparse import csr_matrix
from spatialdata import SpatialData
from spatialdata.models import TableModel

log = logging.getLogger("synthetic_cohort")


# ────────────────────────────────────────────────────────────────────────
# Gene panel + cell-type markers
# ────────────────────────────────────────────────────────────────────────
GENE_PANEL = [
    # Fibrosis axis
    "TGFB1", "TGFBR1", "TGFBR2", "COL1A1", "COL3A1", "ACTA2",
    "PDGFB", "PDGFRB", "FAP",
    # Immune
    "CD80", "CD86", "CTLA4", "CD68", "CD163", "IL6", "IL6R",
    "CXCL12", "CXCR4",
    # Epithelial
    "CDH1", "EGFR", "AGER", "PDPN", "MUC1", "ADAM17",
    # Endothelial
    "VWF", "PLVAP", "LYVE1",
    # Liver-specific
    "ALB", "APOA1", "LRAT", "HNF4A",
    # Housekeeping
    "B2M", "GAPDH", "ACTB", "MALAT1",
    # Filler so panel ~50 genes (named realistically)
    "MKI67", "CDK1", "TOP2A", "FOXJ1", "EPCAM", "KRT19",
    "CD3D", "MS4A1", "MZB1", "CCL19", "CCL21", "CXCL13", "ISG15",
]

CELLTYPE_MARKERS = {
    # Lung
    "AT1":            ["AGER", "PDPN"],
    "Fibroblast":     ["COL1A1", "COL3A1", "PDGFRB", "FAP"],
    "Macrophage_M1":  ["CD68", "CD80", "CD86", "IL6"],
    # Liver
    "Hepatocyte":     ["ALB", "APOA1", "HNF4A"],
    "HSC_activated":  ["ACTA2", "COL1A1", "PDGFRB"],
    "Endothelial":    ["VWF", "PLVAP"],
}

# Condition-specific perturbations (gene → multiplier)
CONDITION_FX = {
    "fibrotic_TBD":     {"TGFB1": 3.0, "COL1A1": 2.5, "ACTA2": 2.0,
                         "ADAM17": 2.0, "FAP": 2.5},
    "fibrotic_nonTBD":  {"TGFB1": 2.0, "COL1A1": 2.0, "ACTA2": 2.0,
                         "FAP": 1.5},
    "control":          {},
}


# ────────────────────────────────────────────────────────────────────────
# Donor manifest for the smoke cohort
# ────────────────────────────────────────────────────────────────────────
SMOKE_DONORS_LUNG = [
    # (subject_id,         tma_core_label, condition,        tbd_status, control_subtype, cores_expected)
    ("SMOKE_LUNG_TBD01",  "TBD1",         "fibrotic_TBD",   "TBD",      None,            2),
    ("SMOKE_LUNG_TBD02",  "TBD2",         "fibrotic_TBD",   "TBD",      None,            2),
    ("SMOKE_LUNG_NS1",    "Control NS1",  "control",        "nonTBD",   "NS",            2),
    ("SMOKE_LUNG_NN1",    "Control NN1",  "control",        "nonTBD",   "NN",            2),
]
SMOKE_DONORS_LIVER = [
    ("SMOKE_LIVER_TBD01", "TBD1",         "fibrotic_TBD",     "TBD",    None,           2),
    ("SMOKE_LIVER_TBD02", "TBD2",         "fibrotic_TBD",     "TBD",    None,           2),
    ("SMOKE_LIVER_CTRL1", "Control 1",    "control",          "nonTBD", "Control",      2),
    ("SMOKE_LIVER_CTRL2", "Control 2",    "control",          "nonTBD", "Control",      2),
    ("SMOKE_LIVER_ALC1",  "Alc Cirh 1",   "fibrotic_nonTBD",  "nonTBD", "AlcCirh",      2),
    ("SMOKE_LIVER_ALC2",  "Alc Cirh 2",   "fibrotic_nonTBD",  "nonTBD", "AlcCirh",      2),
]
LUNG_CELLTYPES  = ["AT1", "Fibroblast", "Macrophage_M1"]
LIVER_CELLTYPES = ["Hepatocyte", "HSC_activated", "Endothelial"]

CELLS_PER_CORE = 250


# ────────────────────────────────────────────────────────────────────────
# Count matrix simulation
# ────────────────────────────────────────────────────────────────────────
def _simulate_core(
    sample_id: str,
    core_id: str,
    tma_core_label: str,
    condition: str,
    celltypes: List[str],
    rng: np.random.Generator,
    n_cells: int = CELLS_PER_CORE,
) -> ad.AnnData:
    """Build one core's AnnData."""
    cells_per_type = np.array_split(np.arange(n_cells), len(celltypes))
    cell_type = np.array(["?"] * n_cells, dtype=object)
    for ct, idx in zip(celltypes, cells_per_type):
        cell_type[idx] = ct

    n_genes = len(GENE_PANEL)
    X = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)

    fx = CONDITION_FX.get(condition, {})

    for i, ct in enumerate(cell_type):
        # Marker genes for this cell type
        for marker in CELLTYPE_MARKERS.get(ct, []):
            j = GENE_PANEL.index(marker)
            X[i, j] += rng.poisson(lam=12.0)
        # Condition perturbation (only on genes affected)
        for gene, mult in fx.items():
            j = GENE_PANEL.index(gene)
            X[i, j] += rng.poisson(lam=2.0 * mult)

    obs = pd.DataFrame({
        "cell_type_L2":   cell_type,
        "core_id":        core_id,
        "tma_core_label": tma_core_label,
    }, index=[f"{sample_id}_cell{i:04d}" for i in range(n_cells)])

    var = pd.DataFrame(index=pd.Index(GENE_PANEL, name="gene"))

    adata = ad.AnnData(
        X=csr_matrix(np.log1p(X)),
        layers={"counts": csr_matrix(X)},
        obs=obs,
        var=var,
    )
    # Random spatial coords in a 500 × 500 µm box
    adata.obsm["spatial"] = rng.uniform(0, 500, size=(n_cells, 2))
    return adata


# ────────────────────────────────────────────────────────────────────────
# TMA builder
# ────────────────────────────────────────────────────────────────────────
def build_tma(
    donors: List[Tuple],
    celltypes: List[str],
    tma_slide: str,
    seed: int = 0,
) -> ad.AnnData:
    """Stack all cores of a TMA into one AnnData."""
    rng = np.random.default_rng(seed)
    parts: List[ad.AnnData] = []
    for sid, label, cond, tbd_st, ctrl_sub, n_cores in donors:
        for core_idx in range(1, n_cores + 1):
            core_id = f"{sid}_core{core_idx:02d}"
            sample_ad = _simulate_core(
                sample_id=core_id,
                core_id=core_id,
                tma_core_label=label,
                condition=cond,
                celltypes=celltypes,
                rng=rng,
            )
            parts.append(sample_ad)

    tma = ad.concat(parts, axis=0, join="outer", merge="same",
                    index_unique=None)
    # TMA-level metadata (sample_id stays unset; concat_samples builds it)
    tma.obs["tma_slide"] = tma_slide
    # Categoricals so downstream tools don't trip on object dtype
    for c in ("core_id", "tma_core_label", "cell_type_L2", "tma_slide"):
        tma.obs[c] = tma.obs[c].astype("category")
    return tma


def write_tma_sdata(adata: ad.AnnData, out_zarr: Path) -> None:
    """Wrap an AnnData into a SpatialData object and write it to zarr."""
    out_zarr.parent.mkdir(parents=True, exist_ok=True)
    if out_zarr.exists():
        import shutil
        shutil.rmtree(out_zarr)
    table = TableModel.parse(adata)
    sdata = SpatialData(tables={"table": table})
    sdata.write(out_zarr)


# ────────────────────────────────────────────────────────────────────────
# Smoke cohort YAML
# ────────────────────────────────────────────────────────────────────────
def write_smoke_cohort_yaml(
    yaml_path: Path,
    lung_tma_sdata: Path,
    liver_tma_sdata: Path,
    project_root: Path,
) -> None:
    """Write a YAML that follows the production schema but points at smoke data.

    Paths to the TMA sdata.zarr are stored RELATIVE TO project_root, matching
    the production yaml convention.
    """
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    relpath = lambda p: str(Path(p).resolve().relative_to(project_root.resolve()))

    donors_yaml = []
    for sid, label, cond, tbd_st, ctrl_sub, n_cores in SMOKE_DONORS_LUNG:
        donors_yaml.append({
            "subject_id": sid,
            "organ": "lung",
            "tma_slide": "LUNG_AF",
            "tma_core_label": label,
            "condition": cond,
            "tbd_status": tbd_st,
            "control_subtype": ctrl_sub,
            "cores_expected": n_cores,
        })
    for sid, label, cond, tbd_st, ctrl_sub, n_cores in SMOKE_DONORS_LIVER:
        donors_yaml.append({
            "subject_id": sid,
            "organ": "liver",
            "tma_slide": "LIVER_AF",
            "tma_core_label": label,
            "condition": cond,
            "tbd_status": tbd_st,
            "control_subtype": ctrl_sub,
            "cores_expected": n_cores,
        })

    payload = {
        "cohort_id": "TBDs_smoke",
        "description": "Synthetic mini-cohort for end-to-end smoke testing of M1+M2",
        "conditions": ["control", "fibrotic_nonTBD", "fibrotic_TBD"],
        "organs": {
            "lung": {
                "panel_name": "Synthetic lung panel (50 genes)",
                "n_genes_total": len(GENE_PANEL),
                "tma_slide": "LUNG_AF",
                "markers_yaml": "pipeline/config/markers/lung_TBDs.yaml",
                "analysis_config": "pipeline/config/config_lung.yaml",
                "tma_raw_dir": "tests/fixtures/_smoke/lung/raw",
                "tma_sdata": relpath(lung_tma_sdata),
                "contrasts": [{"test": "fibrotic_TBD", "ref": "control"}],
            },
            "liver": {
                "panel_name": "Synthetic liver panel (50 genes)",
                "n_genes_total": len(GENE_PANEL),
                "tma_slide": "LIVER_AF",
                "markers_yaml": "pipeline/config/markers/liver_TBDs.yaml",
                "analysis_config": "pipeline/config/config_liver.yaml",
                "tma_raw_dir": "tests/fixtures/_smoke/liver/raw",
                "tma_sdata": relpath(liver_tma_sdata),
                "contrasts": [
                    {"test": "fibrotic_TBD",   "ref": "control"},
                    {"test": "fibrotic_nonTBD", "ref": "control"},
                    {"test": "fibrotic_TBD",   "ref": "fibrotic_nonTBD"},
                ],
            },
        },
        "donors": donors_yaml,
        "integration": {
            "batch_key": "sample_id",
            "n_latent": 10,
            "n_layers": 1,
            "gene_likelihood": "nb",
            "max_epochs": 30,         # tiny for smoke
            "early_stopping": True,
            "use_gpu": False,
        },
        "pseudobulk": {
            "groupby": ["cell_type_L2", "subject_id"],
            "reference_condition": "control",
            "min_cells_per_pseudobulk": 5,
            "min_donors_per_group": 2,
            "fdr_threshold": 0.05,
            "lfc_threshold": 0.585,
        },
        "ccc_multisample": {
            "sample_key": "subject_id",
            "condition_key": "condition",
            "method": "rank_aggregate",
            "resource": "consensus",
            "n_perms": 100,            # fewer perms for speed
            "min_cells": 30,
            "expr_prop": 0.1,
            "fdr_method": "bh",
            "fdr_threshold": 0.05,
        },
        "spacia_meta": {
            "per_sample": True,
            "meta_unit": "subject_id",
            "mcmc_params": "5000,1000,100,2",
            "n_cells_subsample": 200,
            "dist_cutoff_um": 30,
            "multiple_testing": "bonferroni",
            "meta_method": "stouffer",
            "min_donors_validating": 1,
        },
        "pseudotime_cohort": {
            "target_lineages": {
                "lung":  [{"name": "Smoke_lineage", "anchor_celltypes": ["AT1", "Fibroblast"]}],
                "liver": [{"name": "Smoke_lineage", "anchor_celltypes": ["Hepatocyte", "HSC_activated"]}],
            },
            "method": "slingshot",
            "condition_test": "tradeSeq",
        },
        "spatial_domains": {
            "cross_sample_mode": "novae",
            "novae_finetune_epochs": 5,
            "per_core": True,
        },
    }
    with open(yaml_path, "w") as f:
        yaml.safe_dump(payload, f, sort_keys=False)


# ────────────────────────────────────────────────────────────────────────
# One-shot builder
# ────────────────────────────────────────────────────────────────────────
def build_smoke_cohort(out_root: Path, project_root: Path) -> Dict[str, Path]:
    """Create the entire smoke fixture under ``out_root``. Returns key paths."""
    out_root = Path(out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    lung_zarr = out_root / "lung" / "results" / "sdata_TMA.zarr"
    liver_zarr = out_root / "liver" / "results" / "sdata_TMA.zarr"

    log.info("simulating LUNG_AF TMA")
    lung_ad = build_tma(SMOKE_DONORS_LUNG, LUNG_CELLTYPES,
                        tma_slide="LUNG_AF", seed=0)
    write_tma_sdata(lung_ad, lung_zarr)
    log.info("wrote %s (n_obs=%d)", lung_zarr, lung_ad.n_obs)

    log.info("simulating LIVER_AF TMA")
    liver_ad = build_tma(SMOKE_DONORS_LIVER, LIVER_CELLTYPES,
                         tma_slide="LIVER_AF", seed=1)
    write_tma_sdata(liver_ad, liver_zarr)
    log.info("wrote %s (n_obs=%d)", liver_zarr, liver_ad.n_obs)

    yaml_path = out_root / "cohort_smoke.yaml"
    write_smoke_cohort_yaml(yaml_path, lung_zarr, liver_zarr, project_root)
    log.info("wrote %s", yaml_path)

    return {
        "lung_tma": lung_zarr,
        "liver_tma": liver_zarr,
        "cohort_yaml": yaml_path,
        "out_root": out_root,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    repo = Path(__file__).resolve().parents[2]
    out = repo / "tests" / "fixtures" / "_smoke"
    build_smoke_cohort(out, project_root=repo)
    print("\nFixture written to", out)
