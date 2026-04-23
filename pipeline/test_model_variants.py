#!/usr/bin/env python3
"""Test different Cellpose models."""

import numpy as np
from cellpose import models
import traceback

# Create test tile
test_img = np.random.randint(0, 200, (1024, 1024), dtype=np.uint8)
print(f"Test image: shape {test_img.shape}, dtype {test_img.dtype}\n")

# Test different models
for model_name in ['nuclei', 'cyto', 'cyto2', 'cyto3']:
    print(f"Testing model: {model_name:10}", end="", flush=True)
    try:
        m = models.CellposeModel(gpu=True, pretrained_model=model_name)
        print(" ✓ loaded   ", end="", flush=True)

        masks, _, _ = m.eval(test_img, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
        print(f"✓ eval OK - {int(masks.max())} cells")

    except Exception as e:
        error_msg = str(e)[:80]
        print(f"✗ {type(e).__name__}: {error_msg}")
