#!/usr/bin/env python
"""M1 — Pseudobulk DGE via decoupler-py + PyDESeq2.

Per (cell_type, subject_id) the script aggregates counts to a single
pseudobulk profile, then runs PyDESeq2 per cell type with ``condition`` as
the main factor. ``subject_id`` is preserved as the unit of replication
(matches the TMA design: each donor contributes ≥3 technical cores that
already share biology — pooling at donor level is the statistically correct
choice).

Two contrasts (lung + liver), pulled from the cohort YAML under
``organs.<organ>.contrasts``. Each contrast produces:

    TBDs/cohort/results/pseudobulk_decoupler/
        <organ>/<celltype>/<contrast>__decoupler.tsv

Columns: gene, log2FoldChange, lfcSE, stat, pvalue, padj, baseMean.

Run after concat_samples + scvi_integrate (the integrated AnnData has the
cell_type labels finalized).

CLI
---
    python pseudobulk_decoupler.py
        --input  TBDs/cohort/results/cohort_integrated.h5ad
        --cohort pipeline/config/cohort_TBDs.yaml
        [--celltype-col cell_type_L2]
        [--organ lung|liver]
        [--min-cells 10] [--min-donors 2]
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

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import Cohort, load_cohort  # noqa: E402
from utils.paths import cohort_results_root, cohort_yaml_path  # noqa: E402

log = logging.getLogger("pseudobulk_decoupler")


# ────────────────────────────────────────────────────────────────────────
# Pseudobulk aggregation
# ────────────────────────────────────────────────────────────────────────
def build_pseudobulk(
    adata: ad.AnnData,
    celltype_col: str,
    sample_col: str = "subject_id",
    condition_col: str = "condition",
    min_cells: int = 10,
) -> ad.AnnData:
    """Aggregate raw counts per (cell_type, subject_id).

    Returns an AnnData where each obs row is one (cell_type, subject_id)
    pseudobulk, with the original condition/organ/etc. metadata preserved.
    """
    import decoupler as dc  # type: ignore

    if "counts" not in adata.layers:
        raise RuntimeError(
            "pseudobulk needs adata.layers['counts'] — concat_samples preserves it"
        )

    log.info("building pseudobulk: groupby=(%s, %s), min_cells=%d",
             celltype_col, sample_col, min_cells)
    # decoupler 2.x API: dc.pp.pseudobulk (replaces 1.x dc.get_pseudobulk).
    # min_cells filtering is post-hoc via psbulk_n_cells in v2.
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=celltype_col,
        layer="counts",
        mode="sum",
        empty=True,
    )
    # decoupler v2 names the cell-count column ``psbulk_cells`` (v1 used
    # ``psbulk_n_cells``); also keep min_counts protection against vacuous
    # pseudobulks that pass the cell threshold but are still all-zero.
    n_col = next(
        (c for c in ("psbulk_cells", "psbulk_n_cells") if c in pdata.obs.columns),
        None,
    )
    if n_col is not None:
        keep_cells = pdata.obs[n_col] >= min_cells
    else:
        keep_cells = pd.Series(True, index=pdata.obs_names)
    if "psbulk_counts" in pdata.obs.columns:
        keep_counts = pdata.obs["psbulk_counts"] > 0
    else:
        keep_counts = pd.Series(True, index=pdata.obs_names)
    pdata = pdata[keep_cells & keep_counts].copy()
    log.info("pseudobulk after psbulk_cells>=%d filter: n_obs=%d (was %d)",
             min_cells, pdata.n_obs, len(keep_cells))
    log.info("pseudobulk shape: n_obs=%d (cell_type × donor), n_vars=%d",
             pdata.n_obs, pdata.n_vars)
    return pdata


# ────────────────────────────────────────────────────────────────────────
# DE testing (PyDESeq2)
# ────────────────────────────────────────────────────────────────────────
def de_one_celltype(
    pdata: ad.AnnData,
    celltype: str,
    contrast: Dict[str, str],
    celltype_col: str,
    condition_col: str = "condition",
    min_donors: int = 2,
) -> Optional[pd.DataFrame]:
    """One DE test for one cell type, one contrast (test vs ref)."""
    from pydeseq2.dds import DeseqDataSet         # type: ignore
    from pydeseq2.ds import DeseqStats            # type: ignore

    test_lvl = contrast["test"]
    ref_lvl = contrast["ref"]

    mask_ct = (pdata.obs[celltype_col] == celltype).values
    sub = pdata[mask_ct].copy()
    if sub.n_obs == 0:
        return None

    # Subset to the two conditions of interest
    keep = sub.obs[condition_col].isin([test_lvl, ref_lvl])
    sub = sub[keep].copy()

    # Donor count per condition
    by_cond = sub.obs[condition_col].value_counts()
    if by_cond.get(test_lvl, 0) < min_donors or by_cond.get(ref_lvl, 0) < min_donors:
        log.info("[%s] %s vs %s: insufficient donors (%s) — skip",
                 celltype, test_lvl, ref_lvl, dict(by_cond))
        return None

    counts_df = pd.DataFrame(
        sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X),
        index=sub.obs_names,
        columns=sub.var_names,
    ).astype(int)
    metadata = sub.obs[[condition_col]].copy()
    metadata[condition_col] = pd.Categorical(
        metadata[condition_col], categories=[ref_lvl, test_lvl]
    )

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors=condition_col,
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()

    ds = DeseqStats(dds, contrast=[condition_col, test_lvl, ref_lvl], quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    res["gene"] = res.index
    res = res.reset_index(drop=True)
    res["celltype"] = celltype
    res["contrast"] = f"{test_lvl}_vs_{ref_lvl}"
    return res


# ────────────────────────────────────────────────────────────────────────
# Orchestration
# ────────────────────────────────────────────────────────────────────────
def run_organ(
    pdata: ad.AnnData,
    cohort: Cohort,
    organ: str,
    celltype_col: str,
    min_donors: int,
    outdir: Path,
) -> List[Path]:
    organ_pd = pdata[pdata.obs["organ"] == organ].copy()
    contrasts = cohort.organs[organ].get("contrasts", [])
    if not contrasts:
        log.warning("[%s] no contrasts declared in cohort YAML", organ)
        return []

    cell_types = sorted(organ_pd.obs[celltype_col].unique())
    log.info("[%s] %d cell types × %d contrasts", organ, len(cell_types), len(contrasts))

    written: List[Path] = []
    for ct in cell_types:
        for contrast in contrasts:
            res = de_one_celltype(
                organ_pd, ct, contrast,
                celltype_col=celltype_col,
                min_donors=min_donors,
            )
            if res is None or res.empty:
                continue
            organ_dir = outdir / organ / _safe(ct)
            organ_dir.mkdir(parents=True, exist_ok=True)
            cname = f"{contrast['test']}_vs_{contrast['ref']}"
            out = organ_dir / f"{cname}__decoupler.tsv"
            res.to_csv(out, sep="\t", index=False)
            written.append(out)
            log.info("[%s/%s] %s → %s (%d genes)",
                     organ, ct, cname, out.name, len(res))
    return written


def _safe(s: str) -> str:
    return s.replace("/", "_").replace(" ", "_")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="cohort h5ad (post-scVI; or post-concat if running before)")
    p.add_argument("--cohort", type=Path, default=cohort_yaml_path())
    p.add_argument("--outdir", type=Path,
                   default=cohort_results_root() / "pseudobulk_decoupler")
    p.add_argument("--celltype-col", default="cell_type_L2")
    p.add_argument("--organ", choices=("lung", "liver"), default=None)
    p.add_argument("--min-cells", type=int, default=None)
    p.add_argument("--min-donors", type=int, default=None)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    pb_cfg = cohort.pseudobulk
    min_cells = args.min_cells or pb_cfg.get("min_cells_per_pseudobulk", 10)
    min_donors = args.min_donors or pb_cfg.get("min_donors_per_group", 2)

    adata = ad.read_h5ad(args.input)
    log.info("loaded %s (n_obs=%d, n_vars=%d)",
             args.input, adata.n_obs, adata.n_vars)

    pdata = build_pseudobulk(
        adata,
        celltype_col=args.celltype_col,
        sample_col="subject_id",
        min_cells=min_cells,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    organs = [args.organ] if args.organ else list(cohort.organs.keys())
    total: List[Path] = []
    for org in organs:
        total.extend(run_organ(
            pdata, cohort, org,
            celltype_col=args.celltype_col,
            min_donors=min_donors,
            outdir=args.outdir,
        ))
    log.info("done — %d DE tables written under %s", len(total), args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
