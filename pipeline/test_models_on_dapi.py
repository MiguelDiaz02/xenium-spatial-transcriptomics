#!/usr/bin/env python3
"""Test different Cellpose models on actual DAPI image."""

import numpy as np
import spatialdata as sd
from cellpose import models
import xarray as xr

# Load DAPI image
sdata = sd.read_zarr("../human_lung_cancer/results/sdata.zarr")
img_element = sdata["morphology_focus"]

if hasattr(img_element, "ds"):
    scale0 = img_element["scale0"].ds
    img_da = next(iter(scale0.data_vars.values()))
else:
    img_da = img_element["scale0"]

img_uint16 = np.array(img_da.data)
if img_uint16.ndim == 3:
    img_uint16 = img_uint16[0]

# Get center tile with signal
tile_uint16 = img_uint16[16583:17607, 26518:27542]
tile_uint8 = (tile_uint16.astype(np.float32) / tile_uint16.max() * 255).astype(np.uint8)

print(f"Tile uint16: shape {tile_uint16.shape}, min/max: {tile_uint16.min()}/{tile_uint16.max()}")
print(f"Tile uint8:  shape {tile_uint8.shape}, min/max: {tile_uint8.min()}/{tile_uint8.max()}\n")

# Test different models
for model_name in ['nuclei', 'cyto', 'cyto2', 'cyto3']:
    print(f"Testing {model_name:10} with uint16: ", end="", flush=True)
    try:
        m = models.CellposeModel(gpu=True, pretrained_model=model_name)
        masks, _, _ = m.eval(tile_uint16, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0, bsize=512)
        print(f"✓ {int(masks.max())} cells")
    except Exception as e:
        print(f"✗ {type(e).__name__}: {str(e)[:60]}")

    print(f"Testing {model_name:10} with uint8:  ", end="", flush=True)
    try:
        m = models.CellposeModel(gpu=True, pretrained_model=model_name)
        masks, _, _ = m.eval(tile_uint8, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0, bsize=512)
        print(f"✓ {int(masks.max())} cells")
    except Exception as e:
        print(f"✗ {type(e).__name__}: {str(e)[:60]}")
    print()
