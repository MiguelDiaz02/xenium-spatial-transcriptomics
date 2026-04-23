#!/usr/bin/env python3
"""Test different preprocessing approaches for DAPI uint16 image."""

import numpy as np
import spatialdata as sd
from cellpose import models

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
tile = img_uint16[16583:17607, 26518:27542]
print(f"Original tile: shape {tile.shape}, dtype {tile.dtype}, min/max: {tile.min()}/{tile.max()}\n")

# Test different preprocessing approaches
approaches = []

# Approach 1: uint8 via linear scaling
tile1 = (tile.astype(np.float32) / tile.max() * 255).astype(np.uint8)
approaches.append(("uint8 linear scale", tile1))

# Approach 2: uint8 via percentile clipping
p99 = np.percentile(tile, 99)
tile2 = np.clip(tile.astype(np.float32), 0, p99)
tile2 = (tile2 / p99 * 255).astype(np.uint8)
approaches.append(("uint8 percentile clip", tile2))

# Approach 3: uint8 via min-max normalization
tile3 = ((tile - tile.min()) / (tile.max() - tile.min()) * 255).astype(np.uint8)
approaches.append(("uint8 min-max norm", tile3))

# Approach 4: Keep as uint16 but set proper data range
tile4 = tile.astype(np.uint16)
approaches.append(("uint16 original", tile4))

# Approach 5: Normalize uint16 to 0-1 range then to uint8
tile5 = (tile.astype(np.float32) / tile.max() * 255).astype(np.uint8)
approaches.append(("uint8 (float conv)", tile5))

# Test each approach
model = models.CellposeModel(gpu=True, pretrained_model='nuclei')

for name, test_tile in approaches:
    print(f"{name:30} shape {test_tile.shape}, dtype {test_tile.dtype}, min/max {test_tile.min()}/{test_tile.max()}: ", end="", flush=True)
    try:
        masks, _, _ = model.eval(test_tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
        print(f"✓ {int(masks.max())} cells")
    except Exception as e:
        error_msg = str(e).split('\n')[0][:60]
        print(f"✗ {error_msg}")
