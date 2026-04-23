#!/usr/bin/env python3
"""Direct test of segmentation function."""

import sys
sys.path.insert(0, "/home/mdiaz/Documents/Xenium_project/proyecto_demo_xenium/pipeline")

import spatialdata as sd
import numpy as np
from pathlib import Path
import logging

# Setup logger
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
log = logging.getLogger(__name__)

# Import segmentation function
from cellpose import models
from shapely.geometry import Polygon
from skimage.measure import find_contours, regionprops
import geopandas as gpd
from spatialdata.models import ShapesModel

# Load SpatialData
log.info("Loading SpatialData...")
sdata = sd.read_zarr("../human_lung_cancer/results/sdata.zarr")
log.info(f"Loaded SpatialData with elements: {list(sdata.images.keys())}, {list(sdata.labels.keys())}")

# Test parameters
params = {
    "gpu": True,
    "model": "cyto",
    "diameter": None,
    "flow_threshold": 0.4,
    "cellprob_threshold": 0.0,
    "use_he": False,
    "tile_size": 1024,
    "tile_overlap": 128,
}

he_image_path = "../human_lung_cancer/supplemental_files/Xenium_V1_Human_Lung_Cancer_FFPE_he_image.ome.tif"

log.info(f"Starting segmentation with params: {params}")

# Extract image
use_he = params.get("use_he", True) and Path(he_image_path).exists() and "he_image" in sdata
if use_he and "he_image" in sdata:
    log.info("Using H&E image")
    img_element = sdata["he_image"]
else:
    log.info("Using DAPI morphology_focus image")
    img_element = sdata["morphology_focus"]

# Extract numpy array
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

# Manual tiling
tile_size = params.get("tile_size", 2048)
overlap = params.get("tile_overlap", 256)
log.info(f"Manual tiling: tile_size={tile_size}px, overlap={overlap}px")

# Load model
gpu = params.get("gpu", True)
model_type = params.get("model", "nuclei")
model = models.CellposeModel(gpu=gpu, pretrained_model=model_type)

# Create tile grid
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

# Process tiles
diameter = params.get("diameter") or None
flow_thr = params.get("flow_threshold", 0.4)
cellprob = params.get("cellprob_threshold", 0.0)

tile_count = 0
for tile_idx, y_start in enumerate(y_starts):
    for x_start in x_starts:
        tile_count += 1
        y_end = min(y_start + tile_size, h)
        x_end = min(x_start + tile_size, w)
        tile_img = img_np[y_start:y_end, x_start:x_end]

        if tile_count <= 5 or tile_count % 100 == 0:
            log.info(f"  Tile {tile_count}/{n_tiles} [{y_start}:{y_end}, {x_start}:{x_end}] shape {tile_img.shape}")

        try:
            masks_tile, _, _ = model.eval(
                tile_img,
                diameter=diameter,
                channels=None,
                flow_threshold=flow_thr,
                cellprob_threshold=cellprob,
                bsize=512,
            )

            # Reindex
            masks_tile[masks_tile > 0] += cell_id_offset
            cell_id_offset = int(masks_tile.max())

            # Place in combined mask
            mask_region = masks_combined[y_start:y_end, x_start:x_end]
            new_cells = masks_tile > 0
            mask_region[new_cells] = masks_tile[new_cells]

        except Exception as e:
            log.error(f"  ERROR processing tile {tile_count}: {e}")
            import traceback
            traceback.print_exc()
            raise

log.info(f"Tiling complete. Total cells detected: {int(masks_combined.max())}")

# Convert masks → polygons
log.info("Converting masks to polygons...")
polygons = []
cell_ids = []
pixel_size = float(img_element.attrs.get("transform", {}).get("scale", [0.2125, 0.2125])[0])

region_count = 0
for region in regionprops(masks_combined):
    region_count += 1
    if region_count % 1000 == 0:
        log.info(f"  Processing region {region_count}...")

    contours = find_contours(masks_combined == region.label, 0.5)
    if not contours:
        continue
    contour = max(contours, key=len)
    coords_um = contour[:, ::-1] * pixel_size
    if len(coords_um) < 4:
        continue
    poly = Polygon(coords_um)
    if poly.is_valid and not poly.is_empty:
        polygons.append(poly)
        cell_ids.append(str(region.label))

log.info(f"Created {len(polygons)} polygons from {region_count} regions")

# Save
gdf = gpd.GeoDataFrame({"geometry": polygons}, index=cell_ids)
shapes = ShapesModel.parse(gdf)
sdata["cell_boundaries_cellpose"] = shapes
log.info(f"Stored {len(polygons)} Cellpose polygons in sdata")

# Write to zarr
log.info("Writing to zarr...")
sdata.write_element("cell_boundaries_cellpose")
log.info("DONE - Step 02 segmentation completed successfully!")
