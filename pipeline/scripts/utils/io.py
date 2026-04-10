"""
Shared I/O helpers for the xenium_pipeline scripts.

All scripts use a single SpatialData zarr store that is updated in place.
Each step writes only the elements it modifies, then touches a done file.
"""

from __future__ import annotations

import sys
from pathlib import Path

import spatialdata as sd


def read_sdata(zarr_path: str | Path) -> sd.SpatialData:
    """Read SpatialData from zarr store."""
    return sd.read_zarr(str(zarr_path))


def write_element(sdata: sd.SpatialData, element_name: str, zarr_path: str | Path) -> None:
    """
    Write a single element back to an existing zarr store.

    Deletes the element from disk first (if present) to allow clean overwrite,
    then writes the current in-memory version.

    Parameters
    ----------
    sdata
        The SpatialData object (must have the element in memory).
    element_name
        Key of the element to write (e.g., "table", "cell_boundaries").
    zarr_path
        Path to the SpatialData zarr store.
    """
    import zarr

    store = zarr.open(str(zarr_path), mode="r+")
    # Determine the zarr group for this element type
    for group_name in ("tables", "shapes", "images", "points", "labels"):
        if group_name in store and element_name in store[group_name]:
            del store[group_name][element_name]
            break

    # Write element using spatialdata's backing store mechanism
    sdata.write_element(element_name)


def touch_done(done_path: str | Path) -> None:
    """Create a sentinel file marking a pipeline step as complete."""
    Path(done_path).parent.mkdir(parents=True, exist_ok=True)
    Path(done_path).touch()
