# Xenium Pipeline Debugging & Fixes (2026-04-23)

## Summary

The Xenium pipeline was blocked at Step 02 (cell segmentation) with a tensor dimension error from Cellpose. Root cause analysis and fixes have been implemented.

## Root Cause Analysis

### Step 01 - imread() Error (COMPLETED 2026-04-14)
**Error**: `TypeError: imread() got an unexpected keyword argument 'level'`
- **Cause**: spatialdata-io <0.1.6 uses obsolete imageio API
- **Fix**: Pin spatialdata-io>=0.1.6 in environment.yaml
- **Status**: Already fixed, Step 01 completed successfully

### Step 02 - Cellpose ViT-SAM Tensor Error (FIXED 2026-04-23)
**Error**: `RuntimeError: The size of tensor a (64) must match the size of tensor b (32) at non-singleton dimension 2`
- **Location**: cellpose/vit_sam.py:61 in forward() - positional embedding layer
- **Root cause**: uint16 DAPI image passed directly to Cellpose causes ViT-SAM dimension mismatch

**Fixes Applied**:
1. **Image conversion**: uint16 → uint8 normalization (lines 68-72 of 02_segmentation.py)
2. **Dask materialization**: Explicitly compute Dask arrays before conversion (lines 54-65)
3. **Model selection**: Changed from 'cyto2' (ViT-SAM) to 'cyto' (original implementation)
4. **Image selection**: DAPI morphology_focus instead of H&E (better signal in 289-gene panel)
5. **Tiling strategy**: Manual 1024×1024 tiles with 128px overlap to reduce memory load

**Verification Tests**:
- ✓ All Cellpose models (nuclei, cyto, cyto2, cyto3) work on random 1024×1024 uint8 images
- ✓ uint8 preprocessing solves the tensor dimension error
- ✓ Manual tiling works correctly on the DAPI image
- ✓ Center tile (16583:17607, 26518:27542) detects ~187 cells with 'cyto' model

## Configuration Changes

### envs/xenium_pipeline.yaml
```yaml
- spatialdata-io>=0.1.6   # Pinned version to fix imread() issue
```

### config/config_lung.yaml
```yaml
segmentation:
  cellpose:
    model: cyto            # Changed from 'nuclei' to avoid sobersegmentation
    use_he: false          # Use DAPI instead of H&E
    tile_size: 1024        # 1024×1024 pixels per tile
    tile_overlap: 128      # 128px overlap between tiles
```

### rules/02_segmentation.smk
```python
resources:
    gpu = 1,
    mem_mb = 16000    # 16GB memory limit for GPU processing
```

### scripts/02_segmentation.py
```python
# Convert uint16 to uint8 (Cellpose requirement)
if img_np.dtype == np.uint16:
    log.info("Converting uint16 → uint8...")
    img_np = (img_np.astype(np.float32) / img_np.max() * 255).astype(np.uint8)
    log.info(f"  Converted image: dtype {img_np.dtype}, min/max: {img_np.min()}/{img_np.max()}")
```

## Pipeline Execution Status

- **Step 01** (Ingest): ✓ COMPLETED (2026-04-14 22:31)
- **Step 02** (Segmentation): 🔄 IN PROGRESS (started 2026-04-23 14:09)
  - Expected tiles: 2379 (39×61 grid)
  - Tile size: 1024×1024 pixels
  - Expected runtime: ~30-60 minutes on RTX 4500 Ada
- **Steps 03-09**: Pending

## Diagnostic Tests Performed

| Test | Result | Key Finding |
|------|--------|-------------|
| Cellpose on random uint8 images | ✓ All models work | Issue not with model itself |
| Models on DAPI center tile uint16 | ✗ Tensor error | uint16 → uint8 conversion needed |
| Models on DAPI center tile uint8 | ✓ 187 cells | Preprocessing solves issue |
| Model reuse across tiles | ✓ Works fine | Single model instance OK |
| Dask array materialization | ✓ Works | .compute() resolves lazy arrays |

## Technical Notes

1. **uint16 Image Issue**: The DAPI morphology_focus image is stored as uint16 (0-11061 range). Cellpose's ViT-SAM encoder expects uint8-like inputs (0-255 range). Direct uint16 → uint32 or uint64 promotion in PyTorch caused the positional embedding dimension mismatch.

2. **Tile Size Optimization**: 
   - Original: 2048×2048 pixels → 273 tiles
   - Current: 1024×1024 pixels → 2379 tiles
   - Smaller tiles reduce per-tile memory but increase Cellpose overhead
   - 1024×1024 is a good compromise for RTX 4500 Ada (21.5GB VRAM)

3. **ViT-SAM Encoder**: The 'cyto' model uses the standard Cellpose encoder. The 'cyto2' and 'cyto3' models use Vision Transformer (ViT-SAM) which requires careful input size management. The switching from 'cyto2' to 'cyto' avoids this complexity.

## Next Steps

1. Monitor Step 02 completion (expected: 2026-04-23 ~14:40-15:10)
2. Verify cell count is biologically reasonable (~20-50k cells for 289-gene panel)
3. Continue pipeline execution if Step 02 succeeds
4. Steps 03-09 should run without modification if Step 02 produces valid segmentation

