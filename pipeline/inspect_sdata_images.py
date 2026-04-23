#!/usr/bin/env python3
"""Inspect all images in SpatialData."""

import spatialdata as sd
import numpy as np

sdata_path = "../human_lung_cancer/results/sdata.zarr"
print(f"Loading {sdata_path}...")
sdata = sd.read_zarr(sdata_path)

print("\n" + "=" * 70)
print("AVAILABLE IMAGES IN SPATIALDATA")
print("=" * 70)

for key in sdata.images.keys():
    print(f"\n[{key}]", flush=True)
    img = sdata.images[key]
    print(f"  Type: {type(img)}")
    print(f"  Dims: {img.dims}")
    print(f"  Shape: {img.shape}")
    print(f"  Dtype: {img.dtype}")

    # Extract data
    import xarray as xr
    if hasattr(img, "ds"):
        scale0_data = img["scale0"].data
    elif isinstance(img, xr.DataArray):
        scale0_data = img.data
    else:
        scale0_data = img["scale0"].data

    # Convert to numpy if needed
    if hasattr(scale0_data, 'compute'):
        scale0_np = scale0_data.compute()
    else:
        scale0_np = np.array(scale0_data)

    # Get stats per channel
    if scale0_np.ndim == 3:
        for c in range(scale0_np.shape[0]):
            ch_data = scale0_np[c]
            print(f"    Channel {c}: shape {ch_data.shape}, dtype {ch_data.dtype}, min/max/mean: {ch_data.min()}/{ch_data.max()}/{ch_data.mean():.1f}")
    else:
        print(f"    Shape: {scale0_np.shape}, dtype: {scale0_np.dtype}, min/max/mean: {scale0_np.min()}/{scale0_np.max()}/{scale0_np.mean():.1f}")

print("\n" + "=" * 70)
