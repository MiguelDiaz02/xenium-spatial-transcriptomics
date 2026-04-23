"""
Step 09 — Downstream Analyses
================================
Two independent sub-modules (each can be toggled in config):

1. Doublet detection (Scrublet)
   - Simulates synthetic doublets from real cell pairs
   - Assigns doublet score and binary label to each cell
   - Output: adata.obs["doublet_score"], adata.obs["predicted_doublet"]
   - NOTE: Results are annotative only; cells are NOT removed here.
     Removal decision left to the analyst.

2. Pseudobulk differential expression (PyDESeq2)
   - Requires ≥3 biological samples (sample_key column in adata.obs)
   - Uses raw counts from adata.layers["counts"]
   - Compares each cell type vs. all others (one-vs-rest)
   - Output: CSV with DEGs per cell type in results/09_downstream/pseudobulk/
   - DISABLED for single-sample demos (pseudobulk: false in config)
"""


import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def run_doublet_detection(adata):
    """Run Scrublet and add doublet scores to adata.obs."""
    import scrublet as scr
    import numpy as np
    import scipy.sparse as sp

    log.info("Running Scrublet doublet detection ...")
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if sp.issparse(counts):
        counts = counts.toarray()

    scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        n_prin_comps=min(30, adata.n_vars - 1),
    )

    adata.obs["doublet_score"]      = doublet_scores
    adata.obs["predicted_doublet"]  = predicted_doublets

    n_doublets = predicted_doublets.sum()
    pct = n_doublets / len(predicted_doublets) * 100
    log.info(f"Scrublet: {n_doublets:,} predicted doublets ({pct:.1f}%)")
    return adata


def run_pseudobulk_de(adata, group_key: str, sample_key: str, output_dir: Path) -> None:
    """
    Pseudobulk DEA with PyDESeq2.
    Aggregates raw counts per (sample × cell_type), then runs DESeq2.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    output_dir.mkdir(parents=True, exist_ok=True)

    if sample_key not in adata.obs:
        log.warning(f"sample_key '{sample_key}' not in adata.obs — skipping pseudobulk.")
        return

    counts_matrix = adata.layers["counts"]
    if sp.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()

    cell_types = adata.obs[group_key].unique()
    log.info(f"Pseudobulk DE: {len(cell_types)} cell types × samples in '{sample_key}'")

    for ct in cell_types:
        ct_mask = adata.obs[group_key] == ct
        ct_adata = adata[ct_mask]

        # Aggregate raw counts per sample
        samples = ct_adata.obs[sample_key].unique()
        if len(samples) < 3:
            log.warning(f"  {ct}: only {len(samples)} samples — skipping (need ≥3).")
            continue

        agg_counts = {}
        for s in samples:
            s_mask = ct_adata.obs[sample_key] == s
            counts_s = counts_matrix[ct_mask][s_mask].sum(axis=0)
            agg_counts[s] = counts_s

        bulk_df  = pd.DataFrame(agg_counts, index=adata.var_names).T.astype(int)
        metadata = pd.DataFrame({"condition": ["case"] * len(samples)}, index=samples)

        try:
            dds = DeseqDataSet(
                counts=bulk_df,
                metadata=metadata,
                design_factors="condition",
                refit_cooks=True,
            )
            dds.deseq2()
            stat_res = DeseqStats(dds)
            stat_res.summary()
            results_df = stat_res.results_df.sort_values("padj")
            out_path = output_dir / f"de_{ct.replace(' ', '_')}.csv"
            results_df.to_csv(out_path)
            n_sig = (results_df["padj"] < 0.05).sum()
            log.info(f"  {ct}: {n_sig} significant DEGs → {out_path.name}")
        except Exception as exc:
            log.warning(f"  {ct}: DESeq2 failed ({exc})")


def main():
    import spatialdata as sd

    sdata_path  = snakemake.input.sdata                   # noqa: F821
    done_path   = snakemake.output.done                   # noqa: F821
    do_doublets = snakemake.params.doublet_detection      # noqa: F821
    do_pseudo   = snakemake.params.pseudobulk             # noqa: F821
    group_key   = snakemake.params.pseudobulk_group_key   # noqa: F821
    sample_key  = snakemake.params.pseudobulk_sample_key  # noqa: F821

    outdir      = Path(done_path).parent

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.tables[list(sdata.tables.keys())[0]]  # Get first (usually only) table

    if do_doublets:
        adata = run_doublet_detection(adata)
    else:
        log.info("Doublet detection: skipped (doublet_detection: false).")

    if do_pseudo:
        pseudo_dir = outdir / "pseudobulk"
        run_pseudobulk_de(adata, group_key, sample_key, pseudo_dir)
    else:
        log.info("Pseudobulk DE: skipped (pseudobulk: false — requires ≥3 samples).")

    table_name = list(sdata.tables.keys())[0]
    sdata.tables[table_name] = adata
    log.info("Writing downstream results to zarr store ...")
    sdata.delete_element_from_disk(table_name)
    sdata.write_element(table_name)

    Path(done_path).touch()
    log.info("Step 09 — Downstream: DONE")


main()
