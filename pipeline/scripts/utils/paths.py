"""Centralized path resolution for the Xenium pipeline.

Resolution order for the project root:
    1. ``XENIUM_PROJECT_ROOT`` environment variable (explicit override).
    2. Repo layout walk-up: ``<root>/pipeline/scripts/utils/paths.py``.

Dataset selection (e.g. ``human_lung_cancer`` vs ``TBDs_lung_sample01``)
follows the same pattern, controlled by ``XENIUM_DATASET`` (default
``human_lung_cancer`` for backward compatibility with the pilot).

Spacia binary location is resolved from ``SPACIA_PATH`` if set, otherwise
from the in-repo submodule at ``pipeline/external/Spacia/spacia.py``.
"""

from __future__ import annotations

import os
from pathlib import Path

_DEFAULT_DATASET = "human_lung_cancer"


def project_root() -> Path:
    env = os.environ.get("XENIUM_PROJECT_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[3]


def dataset_root(dataset: str | None = None) -> Path:
    name = dataset or os.environ.get("XENIUM_DATASET", _DEFAULT_DATASET)
    return project_root() / name


def results_root(dataset: str | None = None) -> Path:
    return dataset_root(dataset) / "results"


def sdata_path(dataset: str | None = None) -> Path:
    return results_root(dataset) / "sdata.zarr"


def spacia_tool_path() -> Path:
    env = os.environ.get("SPACIA_PATH")
    if env:
        return Path(env).expanduser().resolve()
    return project_root() / "pipeline" / "external" / "Spacia" / "spacia.py"


# ─── P1 multi-sample paths (TBDs cohort) ────────────────────────────────
def cohort_root() -> Path:
    """``$XENIUM_PROJECT_ROOT/TBDs/cohort`` — integrated outputs land here."""
    return project_root() / "TBDs" / "cohort"


def cohort_results_root() -> Path:
    return cohort_root() / "results"


def sample_sdata_path(sample_id: str, organ: str) -> Path:
    """Per-sample SpatialData path: ``TBDs/<organ>/<sample_id>/results/sdata.zarr``."""
    return project_root() / "TBDs" / organ / sample_id / "results" / "sdata.zarr"


def sample_results_dir(sample_id: str, organ: str) -> Path:
    return project_root() / "TBDs" / organ / sample_id / "results"


def cohort_yaml_path() -> Path:
    """Default cohort YAML; overridable per-call via CLI ``--cohort``."""
    return project_root() / "pipeline" / "config" / "cohort_TBDs.yaml"
