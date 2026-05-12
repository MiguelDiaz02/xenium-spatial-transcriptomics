#!/usr/bin/env python
"""M2 — Multi-sample CCC with LIANA+ rank_aggregate + per-condition aggregation.

Extends the pilot F0c_ccc_liana_granular_v2.py (which assumed one sdata.zarr
covering one sample) to the cohort use-case: aggregate LIANA+ results at the
donor level, then contrast conditions to find condition-specific interactions
(e.g., L–R pairs significantly different between fibrotic_TBD and control).

Strategy:
    1. For each donor (subject_id), run LIANA+ rank_aggregate (5-method
       consensus) on the donor's cells.
    2. Stack results into a long table with subject_id, condition.
    3. Aggregate per (cell_type_pair, L_R_pair):
           - n_donors_significant
           - mean magnitude, median rank
           - per-condition mean (test) vs (ref)
    4. Statistical contrast: Mann-Whitney U on each donor's interaction score
       per (cell_type_pair, L_R, condition_pair); BH-FDR across all tests.

Outputs:
    TBDs/cohort/results/ccc_liana_multisample/
        per_donor/<subject_id>__liana.tsv
        aggregated_<organ>.tsv                  — long stacked
        contrast_<organ>__<test>_vs_<ref>.tsv   — Mann-Whitney + FDR
        figures/dotplot_top_<organ>__<contrast>.png

CLI
---
    python F0c_ccc_liana_multisample.py
        --input  TBDs/cohort/results/cohort_integrated.h5ad
        --cohort pipeline/config/cohort_TBDs.yaml
        [--organ lung|liver]
        [--celltype-col cell_type_L2]
        [--n-perms 1000] [--min-cells 50] [--expr-prop 0.1]
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
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import Cohort, load_cohort  # noqa: E402
from utils.paths import cohort_results_root, cohort_yaml_path  # noqa: E402

log = logging.getLogger("liana_multisample")


# ────────────────────────────────────────────────────────────────────────
def _filter_resource_to_panel(resource_df: pd.DataFrame,
                              var_names: pd.Index) -> pd.DataFrame:
    """Keep only L–R pairs whose ligand + receptor are in the gene panel.

    LIANA's consensus resource has ~5k+ pairs covering the whole
    interactome. With sparse panels (Xenium 289/439 genes, our smoke 48
    genes) >98% of the resource is absent, which triggers LIANA's
    'too few features' assertion. Pre-filtering avoids that without
    changing LIANA's defaults for real-data runs.
    """
    panel = set(var_names)
    # LIANA consensus resource columns vary across versions. Defensive lookup:
    lig_col = next(c for c in ("ligand_complex", "ligand") if c in resource_df.columns)
    rec_col = next(c for c in ("receptor_complex", "receptor") if c in resource_df.columns)

    def all_in(col: str) -> pd.Series:
        # Complex columns use '_' to join subunits (e.g. TGFB1_TGFBR1);
        # single-gene columns are atomic — splitting is a no-op in that case.
        return resource_df[col].fillna("").apply(
            lambda s: bool(s) and all(g in panel for g in str(s).split("_"))
        )
    keep = all_in(lig_col) & all_in(rec_col)
    return resource_df.loc[keep].copy()


def run_liana_per_donor(
    adata: ad.AnnData,
    celltype_col: str,
    n_perms: int,
    min_cells: int,
    expr_prop: float,
    resource: str = "consensus",
) -> Optional[pd.DataFrame]:
    """One donor's LIANA+ rank_aggregate result.

    Returns the full long table from LIANA (one row per ordered cell-type pair
    × ligand-receptor pair), or None if too few cells of any cell type or no
    panel-compatible L–R pairs.
    """
    import liana as li  # type: ignore

    # LIANA requires ≥min_cells per cluster; pre-filter
    counts = adata.obs[celltype_col].value_counts()
    eligible = counts[counts >= min_cells].index.tolist()
    if len(eligible) < 2:
        return None
    sub = adata[adata.obs[celltype_col].isin(eligible)].copy()

    # Resource filter — keep only pairs whose every subunit is in the panel.
    resource_df = li.rs.select_resource(resource)
    filtered = _filter_resource_to_panel(resource_df, sub.var_names)
    if filtered.empty:
        log.warning("no L–R pairs in '%s' resource match panel of %d genes",
                    resource, sub.n_vars)
        return None
    log.info("resource '%s' filtered to %d / %d L–R pairs covered by panel",
             resource, len(filtered), len(resource_df))
    # use_raw=False → LIANA reads from sub.X (log1p), matching the convention
    # set by concat_samples (sub.layers['counts'] holds raw counts but LIANA
    # expects log-normalized expression for L–R scoring).
    li.mt.rank_aggregate(
        sub,
        groupby=celltype_col,
        resource=filtered,            # use filtered DataFrame directly
        expr_prop=expr_prop,
        min_cells=min_cells,
        n_perms=n_perms,
        seed=0,
        use_raw=False,
        return_all_lrs=False,
        verbose=False,
    )
    res = sub.uns["liana_res"].copy()
    return res


def aggregate_per_condition(
    long_df: pd.DataFrame, condition_col: str = "condition",
) -> pd.DataFrame:
    """Aggregate per (source, target, ligand_complex, receptor_complex, condition)."""
    grp_cols = ["source", "target", "ligand_complex", "receptor_complex"]
    agg = long_df.groupby(grp_cols + [condition_col]).agg(
        n_donors=("subject_id", "nunique"),
        mean_magnitude=("magnitude_rank", "mean"),
        median_magnitude=("magnitude_rank", "median"),
        mean_specificity=("specificity_rank", "mean"),
    ).reset_index()
    return agg


def contrast_mannwhitney(
    long_df: pd.DataFrame,
    test_lvl: str,
    ref_lvl: str,
    condition_col: str = "condition",
) -> pd.DataFrame:
    """Per (source, target, L-R) test for difference in rank between conditions."""
    grp_cols = ["source", "target", "ligand_complex", "receptor_complex"]
    rows: List[Dict] = []
    for keys, grp in long_df.groupby(grp_cols):
        a = grp.loc[grp[condition_col] == test_lvl, "magnitude_rank"].dropna().values
        b = grp.loc[grp[condition_col] == ref_lvl, "magnitude_rank"].dropna().values
        if len(a) < 2 or len(b) < 2:
            continue
        u, p = st.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({
            **dict(zip(grp_cols, keys)),
            "n_donors_test": len(a),
            "n_donors_ref": len(b),
            "mean_rank_test": float(np.mean(a)),
            "mean_rank_ref": float(np.mean(b)),
            "U": float(u),
            "pvalue": float(p),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["padj"] = multipletests(out["pvalue"], method="fdr_bh")[1]
    return out.sort_values("padj") if not out.empty else out


# ────────────────────────────────────────────────────────────────────────
def run_organ(
    adata: ad.AnnData,
    cohort: Cohort,
    organ: str,
    celltype_col: str,
    outdir: Path,
    n_perms: int,
    min_cells: int,
    expr_prop: float,
) -> None:
    sub = adata[adata.obs["organ"] == organ].copy()
    if sub.n_obs == 0:
        log.warning("[%s] no cells", organ)
        return

    donors = sorted(sub.obs["subject_id"].unique())
    log.info("[%s] running LIANA+ across %d donors", organ, len(donors))

    per_donor_dir = outdir / "per_donor"
    per_donor_dir.mkdir(parents=True, exist_ok=True)

    long_rows: List[pd.DataFrame] = []
    for sid in donors:
        donor_ad = sub[sub.obs["subject_id"] == sid].copy()
        if donor_ad.n_obs < min_cells * 2:
            log.info("[%s] %s: only %d cells — skip", organ, sid, donor_ad.n_obs)
            continue
        res = run_liana_per_donor(
            donor_ad,
            celltype_col=celltype_col,
            n_perms=n_perms,
            min_cells=min_cells,
            expr_prop=expr_prop,
        )
        if res is None or res.empty:
            log.info("[%s] %s: LIANA produced no result", organ, sid)
            continue
        res["subject_id"] = sid
        res["condition"] = str(donor_ad.obs["condition"].iloc[0])
        res["organ"] = organ
        res.to_csv(per_donor_dir / f"{sid}__liana.tsv", sep="\t", index=False)
        long_rows.append(res)
        log.info("[%s] %s: %d LR pairs", organ, sid, len(res))

    if not long_rows:
        log.warning("[%s] no donor results — aborting", organ)
        return

    long_df = pd.concat(long_rows, ignore_index=True)
    long_df.to_csv(outdir / f"aggregated_{organ}.tsv", sep="\t", index=False)

    agg = aggregate_per_condition(long_df)
    agg.to_csv(outdir / f"aggregated_per_condition_{organ}.tsv",
               sep="\t", index=False)

    for cc in cohort.organs[organ].get("contrasts", []):
        contrast = contrast_mannwhitney(long_df, cc["test"], cc["ref"])
        if contrast.empty:
            continue
        out = outdir / f"contrast_{organ}__{cc['test']}_vs_{cc['ref']}.tsv"
        contrast.to_csv(out, sep="\t", index=False)
        log.info("[%s] contrast %s vs %s → %d rows (sig: %d at padj<0.05)",
                 organ, cc["test"], cc["ref"], len(contrast),
                 (contrast["padj"] < 0.05).sum())


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--cohort", type=Path, default=cohort_yaml_path())
    p.add_argument("--outdir", type=Path,
                   default=cohort_results_root() / "ccc_liana_multisample")
    p.add_argument("--organ", choices=("lung", "liver"), default=None)
    p.add_argument("--celltype-col", default="cell_type_L2")
    p.add_argument("--n-perms", type=int, default=None)
    p.add_argument("--min-cells", type=int, default=None)
    p.add_argument("--expr-prop", type=float, default=None)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    cfg = cohort.ccc_multisample
    n_perms = args.n_perms or cfg.get("n_perms", 1000)
    min_cells = args.min_cells or cfg.get("min_cells", 50)
    expr_prop = args.expr_prop or cfg.get("expr_prop", 0.1)

    adata = ad.read_h5ad(args.input)
    log.info("loaded %s (n_obs=%d)", args.input, adata.n_obs)
    args.outdir.mkdir(parents=True, exist_ok=True)

    organs = [args.organ] if args.organ else list(cohort.organs.keys())
    for org in organs:
        run_organ(
            adata, cohort, org,
            celltype_col=args.celltype_col,
            outdir=args.outdir,
            n_perms=n_perms, min_cells=min_cells, expr_prop=expr_prop,
        )

    log.info("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
