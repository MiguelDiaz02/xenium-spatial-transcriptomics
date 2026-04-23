#!/usr/bin/env python3
"""Minimal Cellpose test to diagnose model.eval() issues."""

import sys
import traceback
import numpy as np

print("=" * 60, flush=True)
print("CELLPOSE MINIMAL TEST", flush=True)
print("=" * 60, flush=True)

# Test 1: Import cellpose
try:
    from cellpose import models
    print("✓ Cellpose imported", flush=True)
except Exception as e:
    print(f"✗ Cellpose import failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Test 2: Load model
try:
    print("\nLoading cyto model (GPU)...", flush=True)
    model = models.CellposeModel(gpu=True, pretrained_model='cyto')
    print("✓ Model loaded", flush=True)
except Exception as e:
    print(f"✗ Model load failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Test 3: Simple eval on tiny image
try:
    print("\nTesting eval on 256×256 image...", flush=True)
    test_img = np.random.randint(50, 200, (256, 256), dtype=np.uint8)
    masks, flows, styles = model.eval(test_img, diameter=None, channels=None, flow_threshold=0.4, cellprob_threshold=0.0)
    print(f"✓ eval() OK - detected {int(masks.max())} cells", flush=True)
except Exception as e:
    print(f"✗ eval() failed: {e}", flush=True)
    traceback.print_exc()
    sys.exit(1)

# Test 4: Check tile parameter
try:
    import inspect
    sig = inspect.signature(model.eval)
    has_tile = 'tile' in sig.parameters
    print(f"\neval() has 'tile' parameter: {has_tile}", flush=True)
except Exception as e:
    print(f"Could not inspect eval signature: {e}", flush=True)

print("\n" + "=" * 60, flush=True)
print("ALL TESTS PASSED", flush=True)
print("=" * 60, flush=True)
