#!/usr/bin/env python
"""M3 — Spacia per-core + cohort-level meta-analysis.

Status: SCAFFOLD (full interface, dry-run path executable today; the
``--execute`` path requires the spacia conda env and the per-TMA
sdata.zarr objects to exist).

Architecture
------------
Spacia (Bayesian MIL) operates on one (counts, metadata) pair per run, so
the cohort path is:

    for each donor:
        for each core of that donor:
            export counts + metadata for the (sender, receiver) pair of interest
            invoke spacia via F3_spacia_ccc_validation.run_spacia_job
            collect Pathway_betas.csv

    for each (sender, receiver, pathway):
        meta-analyze across cores of the same donor (technical-rep pooling)
        meta-analyze across donors of the same condition
        compare conditions: Wald-like z-test on combined Stouffer p-values

Outputs (when --execute):
    TBDs/cohort/results/spacia_meta/
        per_core/<sample_id>__<sender>__<receiver>/Pathway_betas.csv
        per_donor_<organ>.tsv       — Stouffer-combined across cores
        per_condition_<organ>.tsv   — Stouffer-combined across donors
        contrast_<organ>__<test>_vs_<ref>.tsv

Meta-analysis: Stouffer's z (default), Fisher's combined p, or
inverse-variance weighted (--meta-method).

Dry-run path (default) — prints the plan and a per-donor JSON manifest of
what jobs *would* run. Useful for QA before launching ~28 × 25 ≈ 700 Spacia
jobs.

CLI
---
    python spacia_meta.py
        --cohort   pipeline/config/cohort_TBDs.yaml
        [--input   TBDs/cohort/results/cohort_integrated.h5ad]
        --pairs    pipeline/config/spacia_pairs_TBDs.tsv     # (sender, receiver) rows
        [--organ lung|liver]
        [--meta-method stouffer|fisher|ivw]
        [--execute]                                          # run Spacia (else dry-run)
        [--dry-run-out plan.json]
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import Cohort, SampleSpec, expand_cores, load_cohort  # noqa: E402
from utils.paths import (  # noqa: E402
    cohort_results_root,
    cohort_yaml_path,
    spacia_tool_path,
)

log = logging.getLogger("spacia_meta")


# ────────────────────────────────────────────────────────────────────────
@dataclass
class SpaciaJob:
    sample_id: str
    subject_id: str
    organ: str
    condition: str
    sender: str
    receiver: str
    outdir: str
    counts_path: str          # filled at execute time
    meta_path: str            # filled at execute time


def plan_jobs(
    cohort: Cohort,
    pairs: pd.DataFrame,
    cohort_h5ad: Path,
    outdir: Path,
    organ_filter: Optional[str] = None,
) -> List[SpaciaJob]:
    """Build the full job list (one job per core × (sender, receiver))."""
    import anndata as ad

    if not cohort_h5ad.exists():
        raise FileNotFoundError(
            f"cohort integrated h5ad missing: {cohort_h5ad} — run "
            f"concat_samples.py first"
        )

    adata = ad.read_h5ad(cohort_h5ad, backed="r")

    # Expand each TMA's donors → cores
    samples: List[SampleSpec] = []
    for org in cohort.organs:
        if organ_filter and org != organ_filter:
            continue
        tma_obs = adata.obs[adata.obs["organ"] == org][["core_id", "tma_core_label"]]
        samples += expand_cores(cohort, tma_obs,
                                tma_slide=cohort.organs[org]["tma_slide"])

    jobs: List[SpaciaJob] = []
    for spec in samples:
        for _, row in pairs.iterrows():
            sender, receiver = row["sender"], row["receiver"]
            jdir = outdir / "per_core" / f"{spec.sample_id}__{sender}__{receiver}"
            jobs.append(SpaciaJob(
                sample_id=spec.sample_id,
                subject_id=spec.subject_id,
                organ=spec.organ,
                condition=spec.condition,
                sender=sender,
                receiver=receiver,
                outdir=str(jdir),
                counts_path="",   # filled at execute time
                meta_path="",
            ))
    return jobs


# ────────────────────────────────────────────────────────────────────────
def stouffer_combine(
    pvals: np.ndarray, weights: Optional[np.ndarray] = None
) -> float:
    """Stouffer's Z combining one-sided p-values."""
    import scipy.stats as st
    if len(pvals) == 0:
        return np.nan
    z = st.norm.isf(pvals)
    if weights is None:
        weights = np.ones_like(z)
    z_comb = np.sum(weights * z) / np.sqrt(np.sum(weights ** 2))
    return float(st.norm.sf(z_comb))


def fisher_combine(pvals: np.ndarray) -> float:
    import scipy.stats as st
    if len(pvals) == 0:
        return np.nan
    return float(st.combine_pvalues(pvals, method="fisher")[1])


def meta_across_cores(
    core_results: List[pd.DataFrame],
    meta_method: str = "stouffer",
) -> pd.DataFrame:
    """Aggregate Pathway_betas.csv rows across cores of a donor.

    Input rows must carry columns: Pathway, Beta, pval (from Spacia).
    Output: one row per pathway with meta_pval + mean_beta + n_cores.
    """
    if not core_results:
        return pd.DataFrame()
    long = pd.concat(core_results, ignore_index=True)
    grouped: List[Dict] = []
    for pathway, grp in long.groupby("Pathway"):
        ps = grp["pval"].dropna().values
        if meta_method == "stouffer":
            mp = stouffer_combine(ps)
        elif meta_method == "fisher":
            mp = fisher_combine(ps)
        elif meta_method == "ivw":
            mp = stouffer_combine(ps)   # IVW falls back to stouffer here
        else:
            raise ValueError(f"unknown meta_method {meta_method!r}")
        grouped.append({
            "Pathway": pathway,
            "n_cores": len(ps),
            "mean_beta": float(grp["Beta"].mean()),
            "meta_pval": mp,
        })
    return pd.DataFrame(grouped)


def meta_across_donors(
    donor_results: List[pd.DataFrame], meta_method: str = "stouffer",
) -> pd.DataFrame:
    """Aggregate donor-level meta results across donors of the same condition."""
    if not donor_results:
        return pd.DataFrame()
    long = pd.concat(donor_results, ignore_index=True)
    rows: List[Dict] = []
    for pathway, grp in long.groupby("Pathway"):
        ps = grp["meta_pval"].dropna().values
        if len(ps) == 0:
            continue
        if meta_method == "stouffer":
            mp = stouffer_combine(ps)
        elif meta_method == "fisher":
            mp = fisher_combine(ps)
        else:
            mp = stouffer_combine(ps)
        rows.append({
            "Pathway": pathway,
            "n_donors": len(ps),
            "mean_beta": float(grp["mean_beta"].mean()),
            "condition_meta_pval": mp,
        })
    return pd.DataFrame(rows)


# ────────────────────────────────────────────────────────────────────────
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", type=Path, default=cohort_yaml_path())
    p.add_argument("--input", type=Path,
                   default=cohort_results_root() / "cohort_integrated.h5ad")
    p.add_argument("--pairs", type=Path, required=True,
                   help="TSV with columns: sender, receiver")
    p.add_argument("--outdir", type=Path,
                   default=cohort_results_root() / "spacia_meta")
    p.add_argument("--organ", choices=("lung", "liver"), default=None)
    p.add_argument("--meta-method", choices=("stouffer", "fisher", "ivw"),
                   default="stouffer")
    p.add_argument("--execute", action="store_true",
                   help="actually run Spacia (default: dry-run only)")
    p.add_argument("--dry-run-out", type=Path, default=None)
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    pairs = pd.read_csv(args.pairs, sep="\t")
    for c in ("sender", "receiver"):
        if c not in pairs.columns:
            raise SystemExit(f"--pairs file missing column {c!r}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    jobs = plan_jobs(cohort, pairs, args.input, args.outdir,
                     organ_filter=args.organ)
    log.info("planned %d Spacia jobs (%d cores × %d pairs)",
             len(jobs),
             len({j.sample_id for j in jobs}),
             len(pairs))

    plan = [asdict(j) for j in jobs]
    plan_out = args.dry_run_out or args.outdir / "plan.json"
    plan_out.write_text(json.dumps(plan, indent=2))
    log.info("wrote plan → %s", plan_out)

    if not args.execute:
        log.info("dry-run only — pass --execute to invoke Spacia")
        return 0

    # ─── EXECUTION PATH ─────────────────────────────────────────────────
    # Implementation outline (left as scaffold for the day samples land):
    #   1. Verify env: spacia conda env must be active (R 4.5.3 + Python 3.8).
    #   2. For each job:
    #        a. Subset the per-TMA AnnData to the core.
    #        b. Use F3_spacia_ccc_validation.export_counts_meta(...) to
    #           write counts/metadata tsvs.
    #        c. Call run_spacia_job(sender, receiver, counts, meta, outdir,
    #                              mcmc_params, n_cells, dist_cutoff).
    #        d. Parse Pathway_betas.csv on success.
    #   3. After all cores done:
    #        - groupby(subject_id) → meta_across_cores → per_donor_<organ>.tsv
    #        - groupby(condition) → meta_across_donors → per_condition_<organ>.tsv
    #        - For each (test, ref) in cohort.organs[organ].contrasts:
    #            compute z-test on Stouffer-combined ps between conditions.
    #   4. Apply Bonferroni across all tested pathways (per spacia_meta.multiple_testing).
    raise SystemExit(
        "--execute path not yet wired. Required to flip: import F3_spacia_ccc_validation, "
        "iterate jobs, run, then call meta_across_cores → meta_across_donors → contrast. "
        "See module docstring for the full algorithm; the dry-run output plan.json is "
        "fully usable as input to a per-job runner script."
    )


if __name__ == "__main__":
    raise SystemExit(main())
