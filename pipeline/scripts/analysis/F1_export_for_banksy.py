#!/usr/bin/env python3
"""
F1 helper — export sdata.zarr counts + spatial coords as simple HDF5/CSV
files that the BANKSY R script can load without reticulate.

Outputs:
  counts.h5     (X dense matrix, var_names, obs_names)
  coords.csv    (cell_id, x, y)
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

import spatialdata as sd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.logging import get_logger  # type: ignore

log = get_logger(__name__)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    log.info(f"Loading {args.sdata}")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    log.info(f"  AnnData: {adata.shape}")

    # ── Counts ────────────────────────────────────────────────────────────
    log.info("Densifying counts (cells x genes)...")
    if hasattr(adata.X, "toarray"):
        X = adata.X.toarray().astype(np.float32)
    else:
        X = np.asarray(adata.X, dtype=np.float32)
    log.info(f"  X: {X.shape}, dtype={X.dtype}")

    counts_h5 = args.outdir / "counts.h5"
    log.info(f"Writing {counts_h5} ...")
    with h5py.File(counts_h5, "w") as f:
        # rhdf5 reads transposed by default; store cells x genes (rows=cells)
        f.create_dataset("X", data=X, compression="gzip", compression_opts=4, chunks=True)
        f.create_dataset("obs_names", data=np.array(adata.obs_names.astype(str), dtype="S"))
        f.create_dataset("var_names", data=np.array(adata.var_names.astype(str), dtype="S"))
    log.info(f"  Saved counts.h5 ({counts_h5.stat().st_size / 1e6:.1f} MB)")

    # ── Coords ────────────────────────────────────────────────────────────
    coords = adata.obsm["spatial"]
    df = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        "x": coords[:, 0],
        "y": coords[:, 1],
    })
    coords_csv = args.outdir / "coords.csv"
    df.to_csv(coords_csv, index=False)
    log.info(f"  Saved coords.csv ({len(df):,} rows)")

    log.info(f"\n✓ Export complete in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
