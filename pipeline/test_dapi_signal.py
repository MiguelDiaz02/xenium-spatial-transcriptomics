#!/usr/bin/env python3
"""Test DAPI signal in different regions."""

import spatialdata as sd
import numpy as np

sdata = sd.read_zarr("../human_lung_cancer/results/sdata.zarr")

# Get morphology_focus image
img = sdata["morphology_focus"]

# Extract to numpy
import xarray as xr
if hasattr(img, "ds"):
    scale0 = img["scale0"].ds
    img_da = next(iter(scale0.data_vars.values()))
elif isinstance(img, xr.DataArray):
    img_da = img
else:
    img_da = img["scale0"]

img_np = np.array(img_da.data)
if img_np.ndim == 3:
    img_np = img_np[0]

h, w = img_np.shape
print(f"Full image shape: ({h}, {w})")
print(f"Image dtype: {img_np.dtype}")
print(f"Full image min/max/mean: {img_np.min()}/{img_np.max()}/{img_np.mean():.1f}\n")

# Test different regions
regions = [
    ("top-left", 0, 1024, 0, 1024),
    ("top-center", 0, 1024, w//2 - 512, w//2 + 512),
    ("top-right", 0, 1024, w - 1024, w),
    ("middle-left", h//2 - 512, h//2 + 512, 0, 1024),
    ("middle-center", h//2 - 512, h//2 + 512, w//2 - 512, w//2 + 512),
    ("middle-right", h//2 - 512, h//2 + 512, w - 1024, w),
    ("bottom-left", h - 1024, h, 0, 1024),
    ("bottom-center", h - 1024, h, w//2 - 512, w//2 + 512),
    ("bottom-right", h - 1024, h, w - 1024, w),
]

for name, y0, y1, x0, x1 in regions:
    tile = img_np[y0:y1, x0:x1]
    print(f"{name:20} [{y0:5}-{y1:5}, {x0:6}-{x1:6}]: shape {tile.shape}, min/max/mean = {tile.min():5}/{tile.max():5}/{tile.mean():8.1f}")
