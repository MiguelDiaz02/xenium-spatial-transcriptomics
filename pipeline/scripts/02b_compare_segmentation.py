"""
Step 02b — Segmentation Comparison Report
==========================================
Compares XOA (10x default) vs the active re-segmentation method.
Generates an interactive HTML report with:
  - Cell count per method
  - Distribution of cell area (µm²)
  - Distribution of transcripts per cell
  - Transcript assignment rate
  - Side-by-side spatial overlay (tile from center of tissue)
"""


import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def compute_seg_metrics(sdata, boundary_key: str, transcript_key: str = "transcripts") -> dict:
    """Compute summary metrics for a given segmentation."""
    import geopandas as gpd

    shapes = sdata[boundary_key]
    n_cells = len(shapes)
    areas   = shapes.geometry.area.values   # µm²

    # Transcript assignment: count how many transcripts fall inside boundaries
    if transcript_key in sdata:
        transcripts_df = sdata[transcript_key].compute()
        pts = gpd.GeoDataFrame(
            transcripts_df,
            geometry=gpd.points_from_xy(transcripts_df["x"], transcripts_df["y"]),
        )
        joined = gpd.sjoin(pts, shapes.reset_index(), how="left", predicate="within")
        n_assigned = joined["index_right"].notna().sum()
        assignment_rate = n_assigned / len(transcripts_df)
        tx_per_cell = joined.groupby("index_right").size()
    else:
        assignment_rate = float("nan")
        tx_per_cell = pd.Series(dtype=float)

    return {
        "n_cells": n_cells,
        "area_median": float(np.median(areas)),
        "area_q25": float(np.percentile(areas, 25)),
        "area_q75": float(np.percentile(areas, 75)),
        "transcript_assignment_rate": assignment_rate,
        "tx_per_cell_median": float(tx_per_cell.median()) if len(tx_per_cell) else float("nan"),
        "areas": areas,
        "tx_per_cell": tx_per_cell.values if len(tx_per_cell) else np.array([]),
    }


def build_report(metrics: dict[str, dict], output_path: str, seg_method: str) -> None:
    """Build an interactive HTML report using plotly."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    methods = list(metrics.keys())
    colors  = {"xoa": "#636EFA", "cellpose": "#EF553B", "baysor": "#00CC96"}

    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=[
            "Cell area distribution (µm²)",
            "Transcripts per cell",
            "Summary metrics",
        ],
        horizontal_spacing=0.08,
    )

    # Panel 1: Cell area violin
    for method, m in metrics.items():
        fig.add_trace(
            go.Violin(
                y=m["areas"],
                name=method.upper(),
                line_color=colors.get(method, "gray"),
                box_visible=True,
                meanline_visible=True,
            ),
            row=1, col=1,
        )

    # Panel 2: Transcripts per cell violin
    for method, m in metrics.items():
        if len(m["tx_per_cell"]) > 0:
            fig.add_trace(
                go.Violin(
                    y=m["tx_per_cell"],
                    name=method.upper(),
                    line_color=colors.get(method, "gray"),
                    box_visible=True,
                    meanline_visible=True,
                    showlegend=False,
                ),
                row=1, col=2,
            )

    # Panel 3: Bar chart — n_cells and assignment rate
    fig.add_trace(
        go.Bar(
            x=[m.upper() for m in methods],
            y=[metrics[m]["n_cells"] for m in methods],
            name="N cells",
            marker_color=[colors.get(m, "gray") for m in methods],
            showlegend=False,
        ),
        row=1, col=3,
    )

    fig.update_layout(
        title_text=(
            f"Segmentation Comparison — Active method: <b>{seg_method.upper()}</b><br>"
            f"<sup>Auto-selected based on panel size and median transcripts/cell</sup>"
        ),
        height=500,
        template="plotly_white",
    )

    # Summary table
    rows = []
    for method, m in metrics.items():
        rows.append({
            "Method": method.upper(),
            "N cells": f"{m['n_cells']:,}",
            "Median area (µm²)": f"{m['area_median']:.1f}",
            "Tx assignment rate": f"{m['transcript_assignment_rate']:.1%}",
            "Median tx/cell": f"{m['tx_per_cell_median']:.1f}",
            "Active": "✓" if method == seg_method else "",
        })
    df_summary = pd.DataFrame(rows)

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Segmentation Comparison</title>
</head>
<body>
  <h2>Xenium Segmentation Comparison</h2>
  {fig.to_html(full_html=False, include_plotlyjs='cdn')}
  <h3>Summary Table</h3>
  {df_summary.to_html(index=False, border=1, justify='center')}
</body>
</html>"""

    Path(output_path).write_text(html)
    log.info(f"Report written to {output_path}")


def main():
    import spatialdata as sd

    sdata_path  = snakemake.input.sdata      # noqa: F821
    report_path = snakemake.output.report    # noqa: F821
    seg_method  = snakemake.params.seg_method  # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)

    # Determine which segmentation keys are present
    boundary_keys = {}
    if "cell_boundaries" in sdata:
        boundary_keys["xoa"] = "cell_boundaries"
    if "cell_boundaries_cellpose" in sdata:
        boundary_keys["cellpose"] = "cell_boundaries_cellpose"
    if "cell_boundaries_baysor" in sdata:
        boundary_keys["baysor"] = "cell_boundaries_baysor"

    if not boundary_keys:
        log.warning("No segmentation boundary elements found in SpatialData. Skipping report.")
        Path(report_path).write_text("<html><body>No segmentation data found.</body></html>")
        return

    log.info(f"Computing metrics for: {list(boundary_keys.keys())}")
    metrics = {}
    for method, key in boundary_keys.items():
        log.info(f"  → {method} ({key}) ...")
        metrics[method] = compute_seg_metrics(sdata, key)

    build_report(metrics, report_path, seg_method)
    log.info("Step 02b — Segmentation comparison: DONE")


main()
