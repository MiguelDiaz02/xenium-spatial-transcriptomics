"""
Step 01 — Ingest
================
1. Extracts the Xenium ZIP bundle to disk.
2. Reads the extracted data with spatialdata-io into a SpatialData object.
3. Registers the H&E image using the affine transform from the alignment CSV.
4. Writes the SpatialData zarr store.

The H&E registration uses the 2×3 affine matrix provided by 10x Genomics
in the alignment CSV (columns: he_x, he_y, xen_x, xen_y for fiducial pairs).
"""


import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.logging import get_logger

log = get_logger(__name__, snakemake.log[0])  # noqa: F821


def extract_zip(zip_path: str, target_dir: str) -> Path:
    """Extract zip bundle; skip if already extracted."""
    target = Path(target_dir)
    sentinel = target / "experiment.xenium"
    if sentinel.exists():
        log.info(f"Already extracted → {target}")
        return target

    target.mkdir(parents=True, exist_ok=True)
    log.info(f"Extracting {zip_path} → {target} ...")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(target)
    log.info("Extraction complete.")
    return target


def build_he_transform(alignment_csv: str):
    """
    Parse the 10x H&E alignment CSV and return a spatialdata Scale2D + Translation.

    The alignment CSV contains a 3×3 affine matrix (row-major) that maps
    H&E pixel coordinates → Xenium µm coordinates.
    Columns: transformation_matrix (single cell, JSON-encoded 3×3 list).
    """
    import json
    import numpy as np
    from spatialdata.transformations import Affine

    import pandas as pd
    df = pd.read_csv(alignment_csv, header=None)

    # 10x format: each row is one row of the 3×3 affine matrix
    # Some versions store as a flat list of 9 values; others as 3×3
    values = df.values.flatten().astype(float)
    if len(values) == 9:
        matrix = values.reshape(3, 3)
    else:
        raise ValueError(f"Unexpected alignment CSV shape: {df.shape}")

    # Build 4×4 homogeneous affine (spatialdata uses 3D transformations)
    affine_4x4 = np.eye(4)
    affine_4x4[:2, :2] = matrix[:2, :2]
    affine_4x4[:2, 3]  = matrix[:2, 2]

    return Affine(affine_4x4, input_axes=("x", "y"), output_axes=("x", "y"))


def main():
    import spatialdata as sd
    import spatialdata_io

    zip_path      = snakemake.input.zip          # noqa: F821
    he_image_path = snakemake.input.he_image     # noqa: F821
    he_align_path = snakemake.input.he_align     # noqa: F821
    extracted_dir = snakemake.params.extracted_dir  # noqa: F821
    sdata_path    = snakemake.output.sdata       # noqa: F821
    done_path     = snakemake.output.done        # noqa: F821
    n_jobs        = snakemake.threads            # noqa: F821

    # 1. Extract ZIP
    xenium_dir = extract_zip(zip_path, extracted_dir)

    # 2. Ingest with spatialdata-io
    log.info("Reading Xenium data with spatialdata-io ...")
    sdata = spatialdata_io.xenium(
        path=str(xenium_dir),
        n_jobs=n_jobs,
    )
    log.info(f"SpatialData loaded: {sdata}")

    # 3. Register H&E image
    if Path(he_image_path).exists() and Path(he_align_path).exists():
        log.info("Registering H&E image ...")
        try:
            from spatialdata_io import xenium_aligned_image
            he_element = xenium_aligned_image(
                image_path=he_image_path,
                alignment_file=he_align_path,
            )
            sdata["he_image"] = he_element
            log.info("H&E image registered.")
        except Exception as exc:
            log.warning(f"H&E registration failed ({exc}); continuing without it.")
    else:
        log.warning("H&E image or alignment file not found; skipping registration.")

    # 4. Write SpatialData zarr store
    log.info(f"Writing SpatialData to {sdata_path} ...")
    sdata.write(sdata_path, overwrite=False)
    log.info("Write complete.")

    Path(done_path).touch()
    log.info("Step 01 — Ingest: DONE")


main()
