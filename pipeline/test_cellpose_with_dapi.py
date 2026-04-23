#!/usr/bin/env python3
"""Test Cellpose on actual DAPI image from SpatialData."""

import sys
sys.path.insert(0, str("/home/mdiaz/Documents/Xenium_project/proyecto_demo_xenium/pipeline"))

import spatialdata as sd
import numpy as np
from cellpose import models
import traceback

print("=" * 60)
print("CELLPOSE DAPI IMAGE TEST")
print("=" * 60)

# Load SpatialData
sdata_path = "../human_lung_cancer/results/sdata.zarr"
print(f"\nLoading SpatialData from {sdata_path}...", flush=True)
try:
    sdata = sd.read_zarr(sdata_path)
    print("✓ SpatialData loaded", flush=True)
except Exception as e:
    print(f"✗ Failed to load: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Extract DAPI image
print("\nExtracting DAPI morphology_focus image...", flush=True)
try:
    img_element = sdata["morphology_focus"]
    print(f"✓ Image element found: {img_element}", flush=True)

    # Extract numpy array at full resolution
    import xarray as xr
    if hasattr(img_element, "ds"):
        scale0 = img_element["scale0"].ds
        img_da = next(iter(scale0.data_vars.values()))
    elif isinstance(img_element, xr.DataArray):
        img_da = img_element
    else:
        img_da = img_element["scale0"]

    img_np = np.array(img_da.data)
    print(f"  Raw shape: {img_np.shape}, dtype: {img_np.dtype}", flush=True)

    # Extract first channel if multi-channel
    if img_np.ndim == 3:
        img_np = img_np[0]
        print(f"  Extracted channel 0, shape now: {img_np.shape}", flush=True)

    print(f"✓ Image extracted: shape {img_np.shape}, dtype {img_np.dtype}", flush=True)
except Exception as e:
    print(f"✗ Image extraction failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Load Cellpose model
print("\nLoading Cellpose model...", flush=True)
try:
    model = models.CellposeModel(gpu=True, pretrained_model='cyto')
    print("✓ Model loaded", flush=True)
except Exception as e:
    print(f"✗ Model load failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Test on first tile
print("\nTesting eval on first 1024×1024 tile...", flush=True)
try:
    tile_img = img_np[0:1024, 0:1024]
    print(f"  Tile shape: {tile_img.shape}, dtype: {tile_img.dtype}, min/max: {tile_img.min()}/{tile_img.max()}", flush=True)

    print("  Running model.eval()...", flush=True)
    masks, flows, styles = model.eval(
        tile_img,
        diameter=None,
        channels=None,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
        bsize=512
    )
    print(f"✓ eval() completed - detected {int(masks.max())} cells", flush=True)
except Exception as e:
    print(f"✗ eval() failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 60)
print("TEST PASSED")
print("=" * 60)
