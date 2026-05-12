#!/usr/bin/env python
"""M1 — scVI batch correction across cohort samples.

Trains an scVI model on the pooled cohort with ``batch_key='sample_id'``
and writes back:

    cohort.obsm['X_scVI']         — 30-D batch-corrected latent
    cohort.layers['scvi_normalized'] — model-decoded normalized expression
    cohort.uns['scvi_history']    — ELBO + reconstruction loss per epoch

The model itself is persisted to ``<outdir>/scvi_model/`` so it can be
reloaded for transfer-learning (e.g. project a new TBDs sample onto an
already-trained cohort latent).

Hyperparameters are taken from the cohort YAML's ``integration:`` block
but every value is CLI-overridable. Defaults match the production
recommendations from scvi-tools 1.4 for Xenium-scale data.

CLI
---
    python scvi_integrate.py --input cohort.h5ad --outdir <dir>
        [--n-latent 30] [--n-layers 2] [--gene-likelihood nb]
        [--max-epochs 400] [--no-gpu] [--cohort <yaml>]
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT / "pipeline" / "scripts"))
from utils.cohort import load_cohort  # noqa: E402
from utils.paths import cohort_results_root, cohort_yaml_path  # noqa: E402

log = logging.getLogger("scvi_integrate")


def integrate(
    adata: ad.AnnData,
    batch_key: str = "sample_id",
    n_latent: int = 30,
    n_layers: int = 2,
    gene_likelihood: str = "nb",
    max_epochs: int = 400,
    early_stopping: bool = True,
    use_gpu: bool = True,
) -> ad.AnnData:
    # Imported lazily so the module is importable without scvi installed
    import scvi  # type: ignore

    if "counts" not in adata.layers:
        raise RuntimeError(
            "scVI requires raw counts in adata.layers['counts'] — concat_samples "
            "preserves it; re-run if missing."
        )

    # On a TMA cohort, batch_key="sample_id" treats each core as its own batch
    # (max correction); "subject_id" pools cores per donor first (gentler).
    # subject_id is added as a categorical covariate so the model can learn
    # donor-level drift without collapsing core-level technical batch effects.
    log.info("setting up AnnData (batch_key=%s, categorical_covariate=subject_id)",
             batch_key)
    categorical_cov = (
        ["subject_id"] if "subject_id" in adata.obs.columns
                          and batch_key != "subject_id" else None
    )
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key,
        categorical_covariate_keys=categorical_cov,
    )

    log.info("training scVI: n_latent=%d, n_layers=%d, likelihood=%s, max_epochs=%d",
             n_latent, n_layers, gene_likelihood, max_epochs)
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood=gene_likelihood,
    )
    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        accelerator="gpu" if use_gpu else "cpu",
    )

    log.info("computing latent + normalized expression")
    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_normalized"] = model.get_normalized_expression(
        return_numpy=True
    )
    adata.uns["scvi_history"] = {
        k: np.asarray(v).tolist()
        for k, v in model.history.items()
    }

    # UMAP on the corrected latent — primary embedding for cohort figures
    log.info("computing neighbors + UMAP on X_scVI")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
    sc.tl.umap(adata, min_dist=0.3)

    return adata, model


def _make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="cohort h5ad from concat_samples.py")
    p.add_argument("--outdir", type=Path, default=cohort_results_root())
    p.add_argument("--cohort", type=Path, default=cohort_yaml_path(),
                   help="cohort YAML (for default hyperparams)")
    p.add_argument("--batch-key", default=None,
                   help="override integration.batch_key from cohort YAML")
    p.add_argument("--n-latent", type=int, default=None)
    p.add_argument("--n-layers", type=int, default=None)
    p.add_argument("--gene-likelihood", choices=("nb", "zinb", "poisson"), default=None)
    p.add_argument("--max-epochs", type=int, default=None)
    p.add_argument("--no-early-stopping", action="store_true")
    p.add_argument("--no-gpu", action="store_true")
    p.add_argument("--log-level", default="INFO")
    return p


def main() -> int:
    args = _make_parser().parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cohort = load_cohort(args.cohort)
    cfg = cohort.integration

    adata = ad.read_h5ad(args.input)
    log.info("loaded %s (n_obs=%d, n_vars=%d)",
             args.input, adata.n_obs, adata.n_vars)

    adata, model = integrate(
        adata,
        batch_key=args.batch_key or cfg.get("batch_key", "sample_id"),
        n_latent=args.n_latent or cfg.get("n_latent", 30),
        n_layers=args.n_layers or cfg.get("n_layers", 2),
        gene_likelihood=args.gene_likelihood or cfg.get("gene_likelihood", "nb"),
        max_epochs=args.max_epochs or cfg.get("max_epochs", 400),
        early_stopping=(not args.no_early_stopping) and cfg.get("early_stopping", True),
        use_gpu=(not args.no_gpu) and cfg.get("use_gpu", True),
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_h5 = args.outdir / "cohort_integrated.h5ad"
    adata.write_h5ad(out_h5, compression="gzip")
    log.info("wrote %s", out_h5)

    model_dir = args.outdir / "scvi_model"
    model.save(str(model_dir), overwrite=True)
    log.info("wrote scVI model to %s", model_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
