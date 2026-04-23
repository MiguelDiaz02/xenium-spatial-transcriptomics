#!/usr/bin/env python3
"""Test Cellpose with uint8 conversion and center tile."""

import spatialdata as sd
import numpy as np
from cellpose import models
import traceback

print("=" * 70)
print("CELLPOSE UINT8 TEST WITH CENTER TILE")
print("=" * 70)

# Load SpatialData
sdata = sd.read_zarr("../human_lung_cancer/results/sdata.zarr")
img_element = sdata["morphology_focus"]

# Extract to numpy
import xarray as xr
if hasattr(img_element, "ds"):
    scale0 = img_element["scale0"].ds
    img_da = next(iter(scale0.data_vars.values()))
elif isinstance(img_element, xr.DataArray):
    img_da = img_element
else:
    img_da = img_element["scale0"]

img_uint16 = np.array(img_da.data)
if img_uint16.ndim == 3:
    img_uint16 = img_uint16[0]

print(f"\nOriginal uint16 image: shape {img_uint16.shape}, dtype {img_uint16.dtype}")
print(f"  min/max: {img_uint16.min()}/{img_uint16.max()}")

# Convert to uint8 (standard for Cellpose)
print("\nConverting uint16 → uint8...")
img_uint8 = (img_uint16.astype(np.float32) / img_uint16.max() * 255).astype(np.uint8)
print(f"  Converted image: dtype {img_uint8.dtype}, min/max: {img_uint8.min()}/{img_uint8.max()}")

# Load model
print("\nLoading Cellpose model...", flush=True)
model = models.CellposeModel(gpu=True, pretrained_model='cyto')
print("✓ Model loaded", flush=True)

# Test 1: Zero-signal tile with uint16
print("\n[Test 1] Zero-signal tile (uint16)...")
try:
    tile = img_uint16[0:1024, 0:1024]
    print(f"  Tile: shape {tile.shape}, min/max: {tile.min()}/{tile.max()}")
    masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
    print(f"  ✓ OK - {int(masks.max())} cells")
except Exception as e:
    print(f"  ✗ Failed: {type(e).__name__}: {str(e)[:100]}")

# Test 2: Signal-rich tile with uint16
print("\n[Test 2] Signal-rich tile (uint16)...")
try:
    tile = img_uint16[16583:17607, 26518:27542]
    print(f"  Tile: shape {tile.shape}, min/max: {tile.min()}/{tile.max()}")
    masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
    print(f"  ✓ OK - {int(masks.max())} cells")
except Exception as e:
    print(f"  ✗ Failed: {type(e).__name__}: {str(e)[:100]}")

# Test 3: Zero-signal tile with uint8
print("\n[Test 3] Zero-signal tile (uint8)...")
try:
    tile = img_uint8[0:1024, 0:1024]
    print(f"  Tile: shape {tile.shape}, min/max: {tile.min()}/{tile.max()}")
    masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
    print(f"  ✓ OK - {int(masks.max())} cells")
except Exception as e:
    print(f"  ✗ Failed: {type(e).__name__}: {str(e)[:100]}")

# Test 4: Signal-rich tile with uint8
print("\n[Test 4] Signal-rich tile (uint8)...")
try:
    tile = img_uint8[16583:17607, 26518:27542]
    print(f"  Tile: shape {tile.shape}, min/max: {tile.min()}/{tile.max()}")
    masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
    print(f"  ✓ OK - {int(masks.max())} cells")
except Exception as e:
    print(f"  ✗ Failed: {type(e).__name__}: {str(e)[:100]}")

print("\n" + "=" * 70)
