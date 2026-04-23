#!/usr/bin/env python3
"""Test whether model reuse causes issues."""

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

# Get multiple tiles
tiles = [
    ("tile 1", img_uint16[0:1024, 0:1024]),
    ("tile 2", img_uint16[1024:2048, 0:1024]),
    ("tile 3", img_uint16[16583:17607, 26518:27542]),  # center tile with signal
]

# Convert to uint8
tiles = [(name, (t.astype(np.float32) / t.max() * 255).astype(np.uint8)) for name, t in tiles]

print("=" * 70)
print("TEST 1: REUSE SINGLE MODEL ACROSS TILES")
print("=" * 70)

model = models.CellposeModel(gpu=True, pretrained_model='nuclei')
print("Model loaded\n")

for name, tile in tiles:
    print(f"{name:20} min/max {tile.min():3}/{tile.max():3}: ", end="", flush=True)
    try:
        masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
        print(f"✓ {int(masks.max())} cells")
    except Exception as e:
        error_msg = str(e).split('\n')[0][:60]
        print(f"✗ {error_msg}")

print("\n" + "=" * 70)
print("TEST 2: RELOAD MODEL FOR EACH TILE")
print("=" * 70)

for name, tile in tiles:
    print(f"{name:20} min/max {tile.min():3}/{tile.max():3}: ", end="", flush=True)
    try:
        model = models.CellposeModel(gpu=True, pretrained_model='nuclei')
        masks, _, _ = model.eval(tile, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
        print(f"✓ {int(masks.max())} cells")
    except Exception as e:
        error_msg = str(e).split('\n')[0][:60]
        print(f"✗ {error_msg}")
