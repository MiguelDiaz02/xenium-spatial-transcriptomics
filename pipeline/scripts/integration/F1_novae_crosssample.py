#!/usr/bin/env python
"""M5 — Cross-sample spatial domain harmonization (Novae multi-slide + Banksy match).

Status: SCAFFOLD. Dry-run path executable today; --execute requires Novae +
GPU and the per-core sdata.zarr (TBDs/<organ>/cores/<sample_id>/results/sdata.zarr).

Why this matters
----------------
Banksy and Novae are computed per-section. With 28 donors × 3–5 cores ≈ 100+
sections, naïve per-section domain inference yields uncomparable domain IDs
across the cohort. Two strategies, both implemented as alternative paths:

    A) Novae multi-slide mode (preferred):
       Novae natively supports multi-slide training. Passing all cores in a
       single model gives a shared embedding space; domain labels are
       directly comparable across cores.

    B) Banksy per-section + post-hoc matching:
       Run Banksy per core, then match domain labels across cores by marker
       similarity (cosine similarity of mean expression per domain).

CLI
---
    python F1_novae_crosssample.py
        --cohort  pipeline/config/cohort_TBDs.yaml
        --organ   lung|liver|all
        --mode    novae|banksy_match
        [--input  TBDs/cohort/results/cohort_integrated.h5ad]
        [--n-domains 7] [--finetune-epochs 100]
        [--execute]

Outputs (when --execute):
    TBDs/cohort/results/spatial_domains/<organ>/
        domain_assignments.tsv            — sample_id, cell_id, domain
        domain_centroids.tsv              — domain mean-expression profile
        figures/domain_map_<sample_id>.png

The output schema is consistent across both modes so downstream
visualization and statistics treat domain labels uniformly.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import Cohort, load_cohort  # noqa: E402
from utils.paths import cohort_results_root, cohort_yaml_path  # noqa: E402

log = logging.getLogger("novae_crosssample")


# ────────────────────────────────────────────────────────────────────────
def plan(cohort: Cohort, organ_filter: Optional[str]) -> Dict[str, List[str]]:
    """Enumerate the cores that would be processed per organ."""
    out: Dict[str, List[str]] = {}
    for org, ocfg in cohort.organs.items():
        if organ_filter and org != organ_filter:
            continue
        donors = [d.subject_id for d in cohort.by_organ(org)]
        out[org] = donors
    return out


# ────────────────────────────────────────────────────────────────────────
def execute_novae(cohort: Cohort, organ: str, outdir: Path,
                  n_domains: int, finetune_epochs: int) -> None:
    """Multi-slide Novae path. Real implementation is sketched below.

    Pipeline:
        1. Load each core's sdata.zarr; build a Novae multi-slide dataset.
        2. Initialize from the pretrained Novae base checkpoint.
        3. Fine-tune for ``finetune_epochs`` on the cohort.
        4. Predict domain labels per cell; write back to per-core sdata.zarr
           under obs['novae_domain'].
        5. Compute domain centroids and per-core domain maps.

    Required env: xenium_pipeline (novae ≥ 0.4.0, torch ≥ 2.2, CUDA).
    """
    # import novae  # type: ignore
    # import spatialdata as sdio
    #
    # tma_zarr = cohort.tma_sdata_for(organ)
    # sdata = sdio.read_zarr(tma_zarr)
    # adata = sdata.tables["table"]
    #
    # # Treat each core as a slide
    # slides = []
    # for core_id, sub in adata.obs.groupby("core_id"):
    #     mask = adata.obs["core_id"].values == core_id
    #     slides.append(adata[mask].copy())
    #
    # model = novae.Novae.from_pretrained("MICS-Lab/novae-human-0")
    # model.fine_tune(slides, num_epochs=finetune_epochs, n_domains=n_domains)
    #
    # assignments = []
    # for slide in slides:
    #     model.predict(slide)
    #     for cell_id, dom in zip(slide.obs_names, slide.obs["novae_domain"]):
    #         assignments.append({
    #             "cell_id": cell_id,
    #             "core_id": slide.obs["core_id"].iloc[0],
    #             "domain": int(dom),
    #         })
    # pd.DataFrame(assignments).to_csv(outdir / "domain_assignments.tsv",
    #                                  sep="\t", index=False)
    raise SystemExit(
        f"[{organ}] novae path SCAFFOLD — uncomment the block above when "
        f"novae+GPU env is verified. Required inputs: "
        f"{cohort.tma_sdata_for(organ)}"
    )


def execute_banksy_match(cohort: Cohort, organ: str, outdir: Path,
                         n_domains: int) -> None:
    """Per-core Banksy → post-hoc domain matching by marker similarity.

    Pipeline:
        1. For each core: run Banksy (via the R-side wrapper F1_export_for_banksy.py
           + Banksy.run.R) at fixed lambda + k.
        2. Compute per-domain centroid in gene space.
        3. Match domains across cores: hungarian assignment on cosine
           distance of centroids, anchored to a reference core (largest by cell count).
        4. Output harmonized labels + a similarity matrix figure.
    """
    raise SystemExit(
        f"[{organ}] banksy_match SCAFFOLD — see docstring. Calls "
        f"F1_export_for_banksy.py per core and joins centroids."
    )


# ────────────────────────────────────────────────────────────────────────
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cohort", type=Path, default=cohort_yaml_path())
    p.add_argument("--organ", choices=("lung", "liver", "all"), default="all")
    p.add_argument("--mode", choices=("novae", "banksy_match"), default="novae")
    p.add_argument("--input", type=Path,
                   default=cohort_results_root() / "cohort_integrated.h5ad")
    p.add_argument("--outdir", type=Path,
                   default=cohort_results_root() / "spatial_domains")
    p.add_argument("--n-domains", type=int, default=7)
    p.add_argument("--finetune-epochs", type=int, default=100)
    p.add_argument("--execute", action="store_true")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    organ_filter = None if args.organ == "all" else args.organ
    p_plan = plan(cohort, organ_filter)
    log.info("planned organs: %s", list(p_plan.keys()))
    for org, donors in p_plan.items():
        log.info("  %s: %d donors → ~%d cores",
                 org, len(donors), sum(d.cores_expected for d in cohort.by_organ(org)))

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "plan.json").write_text(json.dumps(p_plan, indent=2))

    if not args.execute:
        log.info("dry-run only — pass --execute to compute (GPU required for novae)")
        return 0

    for org in p_plan:
        org_out = args.outdir / org
        org_out.mkdir(parents=True, exist_ok=True)
        if args.mode == "novae":
            execute_novae(cohort, org, org_out,
                          n_domains=args.n_domains,
                          finetune_epochs=args.finetune_epochs)
        else:
            execute_banksy_match(cohort, org, org_out, n_domains=args.n_domains)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
