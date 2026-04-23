#!/usr/bin/env python3
"""Debug segmentation with minimal logging to find the issue."""

import sys
import traceback

try:
    print("[1/5] Importing modules...", flush=True)
    import spatialdata as sd
    import numpy as np
    from pathlib import Path
    from cellpose import models
    print("  ✓ Imports successful", flush=True)

    print("\n[2/5] Loading SpatialData...", flush=True)
    sdata = sd.read_zarr("../human_lung_cancer/results/sdata.zarr")
    print("  ✓ SpatialData loaded", flush=True)

    print("\n[3/5] Extracting and converting image...", flush=True)
    img_element = sdata["morphology_focus"]
    import xarray as xr
    if hasattr(img_element, "ds"):
        scale0 = img_element["scale0"].ds
        img_da = next(iter(scale0.data_vars.values()))
    else:
        img_da = img_element["scale0"]

    img_data = img_da.data
    if hasattr(img_data, 'compute'):
        print("  Converting Dask array...", flush=True)
        img_data = img_data.compute()

    img_np = np.array(img_data)
    if img_np.ndim == 3:
        img_np = img_np[0]

    print(f"  Original: shape {img_np.shape}, dtype {img_np.dtype}", flush=True)

    # Convert uint16 to uint8
    if img_np.dtype == np.uint16:
        print("  Converting uint16 → uint8...", flush=True)
        img_np = (img_np.astype(np.float32) / img_np.max() * 255).astype(np.uint8)
        print(f"  ✓ Converted: dtype {img_np.dtype}, min/max {img_np.min()}/{img_np.max()}", flush=True)

    print("\n[4/5] Loading Cellpose model...", flush=True)
    model = models.CellposeModel(gpu=True, pretrained_model='cyto')
    print("  ✓ Model loaded", flush=True)

    print("\n[5/5] Testing first tile...", flush=True)
    tile = img_np[0:1024, 0:1024]
    print(f"  Tile shape: {tile.shape}, dtype: {tile.dtype}, min/max: {tile.min()}/{tile.max()}", flush=True)

    print("  Running model.eval()...", flush=True, end="")
    sys.stdout.flush()

    masks, _, _ = model.eval(
        tile,
        diameter=None,
        channels=None,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
        bsize=512
    )

    print(" ✓", flush=True)
    print(f"  ✓ Detected {int(masks.max())} cells", flush=True)

    print("\n✓✓✓ SUCCESS ✓✓✓", flush=True)

except Exception as e:
    print(f"\n✗ ERROR: {type(e).__name__}: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)
