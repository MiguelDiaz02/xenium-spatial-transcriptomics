#!/usr/bin/env python
"""M1 — Concatenate per-core AnnData slices from TMA SpatialData into a cohort AnnData.

Architecture (TMA-aware):

    LUNG_AF TMA  →  TBDs/lung/results/sdata_TMA.zarr   (one sdata, all cores)
    LIVER_AF TMA →  TBDs/liver/results/sdata_TMA.zarr  (one sdata, all cores)

Both TMAs must be **pre-dearrayed**: the per-TMA AnnData (``sdata.tables['table']``)
must carry ``obs.core_id`` and ``obs.tma_core_label`` columns that match the
labels in ``cohort_TBDs.yaml``. Dearraying is upstream of this script (sopa /
Xenium Explorer); see SENDA_DORADA.md for the chosen workflow.

For each donor in the cohort:
    - find every core_id whose tma_core_label matches the donor
    - subset the TMA AnnData to that core
    - assign sample_id = "<subject_id>_coreNN"
    - inject cohort metadata (organ, condition, tbd_status, subject_id, ...)

Output: ``$XENIUM_PROJECT_ROOT/TBDs/cohort/results/cohort.h5ad``
        ``                                       /cohort_lung.h5ad`` (if --organ lung)

The resulting AnnData carries:
    .X                       → log1p (matches per-TMA table)
    .layers['counts']        → raw counts
    .obs.sample_id           → core-level ID (categorical)
    .obs.core_id             → original TMA core ID (categorical)
    .obs.subject_id          → DONOR — primary statistical replicate unit
    .obs.organ               → {lung, liver}
    .obs.tma_slide           → {LUNG_AF, LIVER_AF}
    .obs.condition           → {control, fibrotic_nonTBD, fibrotic_TBD}
    .obs.tbd_status          → {TBD, nonTBD}
    .obs.control_subtype     → {NS, NN, Control, AlcCirh, ''}
    .obs.cell_type_L1/L2/L3  → carried from per-TMA annotation
    .obsm['spatial_<sid>']   → per-core spatial coords (each in its own frame)

CLI
---
    python concat_samples.py [--cohort PATH] [--outdir PATH]
        [--organ lung|liver]
        [--require-all] [--lung-tma PATH] [--liver-tma PATH]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sdio

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import Cohort, SampleSpec, expand_cores, load_cohort  # noqa: E402
from utils.paths import (  # noqa: E402
    cohort_results_root,
    cohort_yaml_path,
    project_root,
)

log = logging.getLogger("concat_samples")


# ────────────────────────────────────────────────────────────────────────
# Loading
# ────────────────────────────────────────────────────────────────────────
def _read_tma_table(tma_zarr: Path) -> ad.AnnData:
    """Pull the ``table`` AnnData out of a per-TMA SpatialData object."""
    if not tma_zarr.exists():
        raise FileNotFoundError(f"TMA sdata.zarr not found: {tma_zarr}")
    sdata = sdio.read_zarr(tma_zarr)
    if "table" not in sdata.tables:
        raise RuntimeError(f"no 'table' element in {tma_zarr}")
    adata = sdata.tables["table"].copy()

    for col in ("core_id", "tma_core_label"):
        if col not in adata.obs.columns:
            raise RuntimeError(
                f"{tma_zarr}: obs.{col!r} missing — dearray the TMA before "
                f"running concat_samples (sopa.utils.dearray or Xenium Explorer)."
            )

    return adata


def _subset_to_core(
    tma_adata: ad.AnnData,
    spec: SampleSpec,
) -> ad.AnnData:
    mask = (tma_adata.obs["core_id"].astype(str) == spec.core_id).values
    if not mask.any():
        raise RuntimeError(
            f"core_id {spec.core_id!r} produced 0 cells in TMA {spec.tma_slide}"
        )
    sub = tma_adata[mask].copy()

    sub.obs["sample_id"] = spec.sample_id
    sub.obs["subject_id"] = spec.subject_id
    sub.obs["organ"] = spec.organ
    sub.obs["tma_slide"] = spec.tma_slide
    sub.obs["condition"] = spec.condition
    sub.obs["tbd_status"] = spec.tbd_status
    sub.obs["control_subtype"] = spec.control_subtype or ""

    # Unique obs_names across the cohort
    sub.obs_names = [f"{spec.sample_id}-{n}" for n in sub.obs_names]
    return sub


# ────────────────────────────────────────────────────────────────────────
# Cohort concatenation
# ────────────────────────────────────────────────────────────────────────
def concat_cohort(
    cohort: Cohort,
    tma_overrides: Optional[Dict[str, Path]] = None,
    organ: Optional[str] = None,
    require_all: bool = False,
) -> ad.AnnData:
    """Read both (or one) TMAs, expand cores, subset, concat.

    Parameters
    ----------
    tma_overrides
        Optional mapping ``{organ: Path}`` overriding the YAML's ``tma_sdata``.
    organ
        Restrict to a single organ.
    require_all
        Raise if any donor produces 0 cores on its TMA.
    """
    organs_to_process = [organ] if organ else list(cohort.organs.keys())

    loaded: List[ad.AnnData] = []
    spatial_per_sample: Dict[str, np.ndarray] = {}
    missing_donors: List[str] = []

    for org in organs_to_process:
        tma_zarr = (
            tma_overrides.get(org) if tma_overrides else
            cohort.tma_sdata_for(org)
        )
        log.info("[%s] reading TMA %s", org, tma_zarr)
        tma_adata = _read_tma_table(tma_zarr)
        log.info("[%s] %d cells, %d genes, %d unique cores",
                 org, tma_adata.n_obs, tma_adata.n_vars,
                 tma_adata.obs["core_id"].nunique())

        samples = expand_cores(
            cohort, tma_adata.obs,
            tma_slide=cohort.organs[org]["tma_slide"],
        )
        log.info("[%s] expanded to %d core-samples across %d donors",
                 org, len(samples),
                 len({s.subject_id for s in samples}))

        donors_observed: set[str] = set()
        for spec in samples:
            sub = _subset_to_core(tma_adata, spec)
            loaded.append(sub)
            donors_observed.add(spec.subject_id)
            if "spatial" in sub.obsm:
                spatial_per_sample[spec.sample_id] = sub.obsm["spatial"].copy()

        # Sanity: every donor declared for this organ must produce ≥1 core
        for d in cohort.by_organ(org):
            if d.subject_id not in donors_observed:
                missing_donors.append(d.subject_id)
                msg = (f"donor {d.subject_id} (label '{d.tma_core_label}') "
                       f"produced 0 cores on {d.tma_slide}")
                if require_all:
                    raise RuntimeError(msg)
                log.warning("[skip] %s", msg)

    if not loaded:
        raise RuntimeError("no core-samples loaded — check TMA dearraying")

    # Drop the obsm["spatial"] key BEFORE concat — frames differ per core.
    for a in loaded:
        if "spatial" in a.obsm:
            del a.obsm["spatial"]

    log.info("concatenating %d core-AnnData objects (inner gene join)", len(loaded))
    combined = ad.concat(
        loaded,
        axis=0,
        join="inner",
        merge="same",
        label="batch",
        keys=[a.obs["sample_id"].iloc[0] for a in loaded],
        index_unique=None,
    )

    # Restore per-core spatial coords as namespaced obsm entries
    for sid, coords in spatial_per_sample.items():
        mask = (combined.obs["sample_id"] == sid).values
        full = np.full((combined.n_obs, 2), np.nan, dtype=np.float64)
        full[mask] = coords
        combined.obsm[f"spatial_{sid}"] = full

    # Lock factor levels (reference first → for DE)
    combined.obs["condition"] = pd.Categorical(
        combined.obs["condition"],
        categories=cohort.conditions,
        ordered=False,
    )
    for cat in ("sample_id", "subject_id", "organ", "tma_slide",
                "tbd_status", "control_subtype"):
        combined.obs[cat] = combined.obs[cat].astype("category")

    combined.uns["cohort_id"] = cohort.cohort_id
    combined.uns["cohort_conditions"] = list(cohort.conditions)
    combined.uns["cohort_donors"] = [d.subject_id for d in cohort.donors]
    if missing_donors:
        combined.uns["cohort_missing_donors"] = sorted(set(missing_donors))

    return combined


# ────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--cohort", type=Path, default=cohort_yaml_path())
    parser.add_argument("--outdir", type=Path, default=cohort_results_root())
    parser.add_argument("--organ", choices=("lung", "liver"), default=None)
    parser.add_argument("--lung-tma", type=Path, default=None,
                        help="override path to LUNG_AF sdata_TMA.zarr")
    parser.add_argument("--liver-tma", type=Path, default=None,
                        help="override path to LIVER_AF sdata_TMA.zarr")
    parser.add_argument("--require-all", action="store_true",
                        help="fail if any donor produces 0 cores")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    log.info("cohort %s: %d donors across %d organs",
             cohort.cohort_id, cohort.n_donors(), len(cohort.organs))

    overrides: Dict[str, Path] = {}
    if args.lung_tma:
        overrides["lung"] = args.lung_tma
    if args.liver_tma:
        overrides["liver"] = args.liver_tma

    combined = concat_cohort(
        cohort,
        tma_overrides=overrides or None,
        organ=args.organ,
        require_all=args.require_all,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.organ}" if args.organ else ""
    out = args.outdir / f"cohort{suffix}.h5ad"
    combined.write_h5ad(out, compression="gzip")
    log.info("wrote %s (n_obs=%d, n_vars=%d, n_samples=%d, n_donors=%d)",
             out, combined.n_obs, combined.n_vars,
             combined.obs["sample_id"].nunique(),
             combined.obs["subject_id"].nunique())

    # JSON sidecar for R consumers
    json_out = args.outdir / "cohort_meta.json"
    cohort.to_json(json_out)
    log.info("wrote %s (for R-side scripts)", json_out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
