"""
Step 02 — Cell Re-segmentation
================================
Runs the auto-selected segmentation method (Cellpose or Baysor) and
replaces the cell/nucleus boundary shapes in the SpatialData store.

Method is resolved at DAG-build time in the Snakefile based on:
  - Cellpose: median tx/cell < 100 OR panel genes < 500
  - Baysor:   median tx/cell >= 100 AND panel genes >= 500
"""


import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


# ─── Cellpose ────────────────────────────────────────────────────────────────

def run_cellpose(sdata, params: dict, he_image_path: str):
    """
    Segment cells using Cellpose with manual tiling for large images.
    Cellpose 4.x doesn't support tile= parameter, so we tile manually.
    Returns updated SpatialData with 'cell_boundaries_cellpose' shape element.
    """
    from cellpose import models
    import spatialdata as sd
    from spatialdata.models import ShapesModel
    from shapely.geometry import Polygon
    import geopandas as gpd

    use_he   = params.get("use_he", True) and Path(he_image_path).exists() and "he_image" in sdata
    gpu      = params.get("gpu", True)
    model_type = params.get("model", "nuclei")
    diameter   = params.get("diameter") or None
    flow_thr   = params.get("flow_threshold", 0.4)
    cellprob   = params.get("cellprob_threshold", 0.0)

    # Select image
    if use_he and "he_image" in sdata:
        log.info("Cellpose: using H&E image")
        img_element = sdata["he_image"]
    else:
        log.info("Cellpose: using DAPI morphology_focus image")
        img_element = sdata["morphology_focus"]

    # Extract numpy array at full resolution (s0).
    import xarray as xr
    if hasattr(img_element, "ds"):
        scale0 = img_element["scale0"].ds
        img_da = next(iter(scale0.data_vars.values()))
    elif isinstance(img_element, xr.DataArray):
        img_da = img_element
    else:
        img_da = img_element["scale0"]
    img_np = np.array(img_da.data)
    if img_np.ndim == 3:
        img_np = img_np[0]

    log.info(f"Image shape: {img_np.shape}, dtype: {img_np.dtype}")

    # Manual tiling: Cellpose 4.x doesn't support tile=True
    tile_size = params.get("tile_size", 2048)
    overlap = params.get("tile_overlap", 256)
    log.info(f"Manual tiling: tile_size={tile_size}px, overlap={overlap}px")

    model = models.CellposeModel(gpu=gpu, pretrained_model=model_type)

    # Create tile grid with overlap
    h, w = img_np.shape
    masks_combined = np.zeros((h, w), dtype=np.uint32)
    cell_id_offset = 0

    # Compute tile boundaries
    y_starts = list(range(0, h, tile_size - overlap))
    if y_starts[-1] + tile_size < h:
        y_starts.append(h - tile_size)
    x_starts = list(range(0, w, tile_size - overlap))
    if x_starts[-1] + tile_size < w:
        x_starts.append(w - tile_size)

    n_tiles = len(y_starts) * len(x_starts)
    log.info(f"Processing {n_tiles} tiles ({len(y_starts)}×{len(x_starts)})")

    # Process each tile
    for tile_idx, y_start in enumerate(y_starts):
        for x_start in x_starts:
            y_end = min(y_start + tile_size, h)
            x_end = min(x_start + tile_size, w)
            tile_img = img_np[y_start:y_end, x_start:x_end]

            log.info(f"  Tile [{y_start}:{y_end}, {x_start}:{x_end}] shape {tile_img.shape}")

            # Run Cellpose on tile
            masks_tile, _, _ = model.eval(
                tile_img,
                diameter=diameter,
                channels=None,
                flow_threshold=flow_thr,
                cellprob_threshold=cellprob,
                bsize=512,  # internal tiling for Cellpose
            )

            # Reindex mask labels to avoid collisions
            masks_tile[masks_tile > 0] += cell_id_offset
            cell_id_offset = int(masks_tile.max())

            # Place in combined mask (overlapping regions: keep largest)
            mask_region = masks_combined[y_start:y_end, x_start:x_end]
            # In overlap zones, keep cell with larger label (from newer tiles)
            new_cells = masks_tile > 0
            mask_region[new_cells] = masks_tile[new_cells]

    masks = masks_combined
    n_cells = int(masks.max())
    log.info(f"Cellpose detected {n_cells:,} cells (after tiling).")

    # Convert masks → polygons for ShapesModel
    # Uses regionprops to get contour for each mask label
    from skimage.measure import find_contours, label as sk_label
    from skimage.measure import regionprops

    polygons = []
    cell_ids = []
    pixel_size = float(img_element.attrs.get("transform", {}).get("scale", [0.2125, 0.2125])[0])

    for region in regionprops(masks):
        contours = find_contours(masks == region.label, 0.5)
        if not contours:
            continue
        contour = max(contours, key=len)
        # contour is (row, col) → convert to (x µm, y µm)
        coords_um = contour[:, ::-1] * pixel_size   # (col, row) * pixel_size
        if len(coords_um) < 4:
            continue
        poly = Polygon(coords_um)
        if poly.is_valid and not poly.is_empty:
            polygons.append(poly)
            cell_ids.append(str(region.label))

    gdf = gpd.GeoDataFrame({"geometry": polygons}, index=cell_ids)
    shapes = ShapesModel.parse(gdf)
    sdata["cell_boundaries_cellpose"] = shapes
    log.info(f"Stored {len(polygons):,} Cellpose polygons in sdata['cell_boundaries_cellpose']")
    return sdata


# ─── Baysor ──────────────────────────────────────────────────────────────────

def run_baysor(sdata, params: dict):
    """
    Segment cells using Baysor (transcript-based).
    Requires Julia + Baysor installed (see envs/xenium_pipeline.yaml).
    Uses XOA nucleus boundaries as spatial prior.
    """
    import subprocess
    import tempfile
    import pandas as pd
    import geopandas as gpd
    from spatialdata.models import ShapesModel
    from shapely.geometry import Polygon

    log.info("Baysor segmentation — writing transcripts to temp CSV ...")

    # Extract transcripts as DataFrame
    transcripts = sdata["transcripts"].compute()   # dask → pandas
    transcripts_csv = Path(tempfile.mkdtemp()) / "transcripts.csv"
    transcripts[["x", "y", "z", "gene"]].to_csv(transcripts_csv, index=False)

    # Export nucleus boundaries as CSV for prior
    prior_csv = transcripts_csv.parent / "nucleus_prior.csv"
    if "nucleus_boundaries" in sdata:
        nuc_df = sdata["nucleus_boundaries"].reset_index()
        nuc_df["cell"] = nuc_df.index
        nuc_df[["cell", "geometry"]].to_csv(prior_csv, index=False)
        prior_arg = ["--prior-segmentation-counts", str(prior_csv)]
    else:
        prior_arg = []

    output_dir = transcripts_csv.parent / "baysor_out"
    output_dir.mkdir(exist_ok=True)

    cmd = [
        "baysor", "run",
        "--x", "x", "--y", "y", "--z", "z", "--gene", "gene",
        "--min-molecules-per-cell", str(params["min_molecules_per_cell"]),
        "--scale", str(params["scale"]),
        "--n-clusters", str(params["n_clusters"]),
        "--prior-segmentation-confidence", str(params["prior_segmentation_confidence"]),
        "-o", str(output_dir),
        str(transcripts_csv),
    ] + prior_arg

    log.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Baysor failed:\n{result.stderr}")

    # Parse Baysor output segmentation CSV
    seg_csv = output_dir / "segmentation.csv"
    seg_df  = pd.read_csv(seg_csv)

    # Build polygons from cell assignment
    polygons = []
    cell_ids = []
    for cell_id, group in seg_df.groupby("cell"):
        if cell_id == 0:
            continue   # unassigned transcripts
        pts = group[["x", "y"]].values
        if len(pts) < 4:
            continue
        from scipy.spatial import ConvexHull
        try:
            hull = ConvexHull(pts)
            poly = Polygon(pts[hull.vertices])
            if poly.is_valid and not poly.is_empty:
                polygons.append(poly)
                cell_ids.append(str(cell_id))
        except Exception:
            continue

    gdf    = gpd.GeoDataFrame({"geometry": polygons}, index=cell_ids)
    shapes = ShapesModel.parse(gdf)
    sdata["cell_boundaries_baysor"] = shapes
    log.info(f"Stored {len(polygons):,} Baysor polygons in sdata['cell_boundaries_baysor']")
    return sdata


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    import spatialdata as sd

    sdata_path   = snakemake.input.sdata       # noqa: F821
    done_path    = snakemake.output.done       # noqa: F821
    method       = snakemake.params.method     # noqa: F821
    cellpose_cfg = snakemake.params.cellpose   # noqa: F821
    baysor_cfg   = snakemake.params.baysor     # noqa: F821
    he_image_path = snakemake.params.he_image  # noqa: F821

    log.info(f"Loading SpatialData from {sdata_path} ...")
    sdata = sd.read_zarr(sdata_path)

    if method == "cellpose":
        sdata = run_cellpose(sdata, cellpose_cfg, he_image_path)
        new_element = "cell_boundaries_cellpose"
    elif method == "baysor":
        sdata = run_baysor(sdata, baysor_cfg)
        new_element = "cell_boundaries_baysor"
    elif method == "xoa":
        log.info("method=xoa: keeping original 10x segmentation, no changes.")
        Path(done_path).touch()
        return
    else:
        raise ValueError(f"Unknown segmentation method: {method!r}")

    # Write only the new shapes element back to the zarr store
    log.info(f"Writing {new_element} to zarr store ...")
    sdata.write_element(new_element)

    Path(done_path).touch()
    log.info(f"Step 02 — Segmentation ({method.upper()}): DONE")


main()
