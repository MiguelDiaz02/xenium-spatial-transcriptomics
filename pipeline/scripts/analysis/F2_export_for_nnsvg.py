#!/usr/bin/env python3
"""
Helper: export counts matrix + spatial coords from SpatialData to Matrix Market files
for nnSVG R script consumption. Subsamples to n_cells_subsample (default 10k) because
nnSVG is a Gaussian process method that doesn't scale to 268k cells.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import spatialdata as sd
import scipy.io
import scipy.sparse as sp


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sdata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--n_cells_subsample", type=int, default=10000,
                   help="Max cells to use (nnSVG is O(n^2); default 10000)")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading sdata from {args.sdata}...")
    sdata = sd.read_zarr(str(args.sdata))
    table_name = list(sdata.tables.keys())[0]
    adata = sdata.tables[table_name]
    print(f"  Full dataset: {adata.n_obs:,} cells x {adata.n_vars} genes")

    # Subsample for nnSVG (Gaussian process methods don't scale to 268k cells)
    if adata.n_obs > args.n_cells_subsample:
        np.random.seed(args.seed)
        idx = np.sort(np.random.choice(adata.n_obs, args.n_cells_subsample, replace=False))
        adata = adata[idx].copy()
        print(f"  Subsampled to {adata.n_obs:,} cells (seed={args.seed})")

    # Use raw counts if available
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X

    # Save as Matrix Market (genes x cells) — natively readable by R Matrix::readMM()
    if sp.issparse(X):
        X_gc = X.T.tocsc()
    else:
        X_gc = sp.csc_matrix(X.T)

    scipy.io.mmwrite(str(args.outdir / "counts_matrix.mtx"), X_gc)
    pd.Series(adata.var_names.tolist()).to_csv(
        args.outdir / "gene_names.csv", index=False, header=False)
    pd.Series(adata.obs_names.tolist()).to_csv(
        args.outdir / "cell_names.csv", index=False, header=False)
    pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"]).to_csv(
        args.outdir / "spatial_coords.csv", index=False)

    print(f"  Saved: counts ({X_gc.shape}), coords ({adata.obsm['spatial'].shape})")
    print("Done.")


if __name__ == "__main__":
    main()
