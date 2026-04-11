"""
Step 03 — Quality Control
==========================
Computes per-cell QC metrics and filters low-quality cells.

Metrics computed:
  - total_counts: total transcripts per cell
  - n_genes_by_counts: number of unique genes detected
  - pct_neg_ctrl: fraction of transcripts from negative control probes
  - cell_area: cell polygon area in µm²

Filtering thresholds are read from config (qc section).
Cells failing any threshold are removed from sdata.table.
The filtered table is written back to the zarr store.
"""


import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def compute_qc_metrics(adata, sdata):
    """Add QC columns to adata.obs."""
    import scanpy as sc
    import numpy as np

    # Standard scanpy QC (total counts, n genes, % mito — no mito in Xenium)
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)

    # Negative control transcripts (genes starting with "BLANK_" or "NegControl")
    neg_ctrl_genes = [g for g in adata.var_names
                      if g.startswith("BLANK_") or g.startswith("NegControl")]
    if neg_ctrl_genes:
        adata.obs["pct_neg_ctrl"] = (
            np.asarray(adata[:, neg_ctrl_genes].X.sum(axis=1)).flatten()
            / np.asarray(adata.X.sum(axis=1)).flatten().clip(1)
            * 100
        )
        log.info(f"Found {len(neg_ctrl_genes)} negative control genes.")
    else:
        adata.obs["pct_neg_ctrl"] = 0.0
        log.info("No negative control genes found.")

    # Cell area from shapes element (if available)
    seg_key = None
    for key in ("cell_boundaries_cellpose", "cell_boundaries_baysor", "cell_boundaries"):
        if key in sdata:
            seg_key = key
            break

    if seg_key:
        shapes = sdata[seg_key]
        area_series = shapes.geometry.area
        # Match by cell index (adata.obs_names should match shape index)
        area_dict = area_series.to_dict()
        adata.obs["cell_area_um2"] = adata.obs_names.map(area_dict).fillna(0).values
        log.info(f"Cell area computed from '{seg_key}'.")
    else:
        adata.obs["cell_area_um2"] = 0.0

    return adata


def filter_cells(adata, params: dict) -> tuple:
    """Apply QC thresholds. Returns (filtered_adata, n_removed)."""
    n_before = adata.n_obs
    mask = (
        (adata.obs["total_counts"]   >= params["min_transcripts"]) &
        (adata.obs["total_counts"]   <= params["max_transcripts"]) &
        (adata.obs["n_genes_by_counts"] >= params["min_genes"])    &
        (adata.obs["pct_neg_ctrl"]   <= params["max_pct_neg_ctrl"])
    )
    adata = adata[mask].copy()
    n_after   = adata.n_obs
    n_removed = n_before - n_after
    log.info(
        f"QC filter: {n_before:,} → {n_after:,} cells "
        f"({n_removed:,} removed, {n_removed/n_before:.1%})"
    )
    return adata, n_removed


def build_qc_report(adata_before, adata_after, report_path: str, params: dict) -> None:
    """Generate QC HTML report with metric distributions."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import pandas as pd

    metrics = [
        ("total_counts",       "Transcripts per cell",     params["min_transcripts"], params["max_transcripts"]),
        ("n_genes_by_counts",  "Genes per cell",           params["min_genes"],       None),
        ("pct_neg_ctrl",       "Negative control (%)",     None,                      params["max_pct_neg_ctrl"]),
        ("cell_area_um2",      "Cell area (µm²)",          None,                      None),
    ]

    fig = make_subplots(rows=2, cols=2, subplot_titles=[m[1] for m in metrics])
    positions = [(1,1),(1,2),(2,1),(2,2)]

    for (col, label, vmin, vmax), (row, col_pos) in zip(metrics, positions):
        if col not in adata_before.obs:
            continue
        vals_before = adata_before.obs[col].values
        vals_after  = adata_after.obs[col].values

        fig.add_trace(
            go.Histogram(x=vals_before, name="before", opacity=0.5,
                         marker_color="#636EFA", showlegend=(row==1 and col_pos==1)),
            row=row, col=col_pos,
        )
        fig.add_trace(
            go.Histogram(x=vals_after, name="after QC", opacity=0.7,
                         marker_color="#EF553B", showlegend=(row==1 and col_pos==1)),
            row=row, col=col_pos,
        )
        if vmin is not None:
            fig.add_vline(x=vmin, line_dash="dash", line_color="green",
                          row=row, col=col_pos)
        if vmax is not None:
            fig.add_vline(x=vmax, line_dash="dash", line_color="red",
                          row=row, col=col_pos)

    fig.update_layout(
        title_text=(
            f"QC Report — {adata_before.n_obs:,} → {adata_after.n_obs:,} cells "
            f"({adata_before.n_obs - adata_after.n_obs:,} removed)"
        ),
        height=700,
        template="plotly_white",
        barmode="overlay",
    )

    Path(report_path).write_text(
        f"<!DOCTYPE html><html><head><meta charset='utf-8'></head><body>"
        f"<h2>Xenium QC Report</h2>"
        f"{fig.to_html(full_html=False, include_plotlyjs='cdn')}"
        f"</body></html>"
    )


def main():
    import spatialdata as sd

    sdata_path  = snakemake.input.sdata    # noqa: F821
    done_path   = snakemake.output.done    # noqa: F821
    report_path = snakemake.output.report  # noqa: F821
    params      = {
        "min_transcripts":  snakemake.params.min_transcripts,   # noqa: F821
        "max_transcripts":  snakemake.params.max_transcripts,   # noqa: F821
        "min_genes":        snakemake.params.min_genes,         # noqa: F821
        "max_pct_neg_ctrl": snakemake.params.max_pct_neg_ctrl,  # noqa: F821
    }

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)
    adata = sdata.table

    log.info("Computing QC metrics ...")
    adata = compute_qc_metrics(adata, sdata)

    adata_before = adata.copy()
    adata_after, _ = filter_cells(adata, params)

    # Build report before updating sdata
    log.info("Building QC report ...")
    build_qc_report(adata_before, adata_after, report_path, params)

    # Update sdata.table with filtered cells and write back
    sdata.table = adata_after
    log.info("Writing filtered table back to zarr store ...")
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")

    Path(done_path).touch()
    log.info("Step 03 — QC: DONE")


main()
