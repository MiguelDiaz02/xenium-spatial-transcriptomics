"""Cohort loader for multi-sample Xenium TMA analyses.

The TBDs cohort is built on TWO Xenium TMA slides (LUNG_AF, LIVER_AF). Each
TMA hosts multiple donors, with several technical-replicate cores per donor.
The cohort YAML enumerates DONORS; cores are resolved at runtime by reading
``obs.core_id`` from the per-TMA SpatialData (after dearraying).

Granularity contract:
    SAMPLE   = one TMA core              (statistical *unit of observation*)
    DONOR    = subject_id                (statistical *unit of replication*)
    SLIDE    = one TMA (one Xenium run)  (technical-batch grouping if needed)

Two main entry points:
    load_cohort(yaml_path) -> Cohort
        Parses + validates the YAML. Returns a Cohort object whose
        ``donors`` list reflects what the YAML declares.

    expand_cores(cohort, tma_sdata_obs) -> List[SampleSpec]
        Given a Cohort + the obs DataFrame of a per-TMA sdata.zarr (already
        dearrayed, with columns ``core_id`` and ``tma_core_label``), expands
        each donor into one SampleSpec per core present on the slide.

Schema is validated by ``tests/test_p1_multisample.py::TestCohortYaml``.
"""
from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd
import yaml

from utils.paths import project_root


_SUBJECT_ID_RE = re.compile(r"^[A-Za-z0-9_]+$")
_VALID_TBD = {"TBD", "nonTBD"}
DEFAULT_COHORT_YAML = "pipeline/config/cohort_TBDs.yaml"


# ────────────────────────────────────────────────────────────────────────
# Dataclasses
# ────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class DonorSpec:
    """A donor row from the YAML — pre-expansion."""

    subject_id: str
    organ: str
    tma_slide: str
    tma_core_label: str
    condition: str
    tbd_status: str
    control_subtype: Optional[str] = None
    cores_expected: int = 0
    notes: str = ""


@dataclass(frozen=True)
class SampleSpec:
    """One TMA core — the unit consumed by concat_samples and downstream.

    ``sample_id`` is built deterministically as ``<subject_id>_core<NN>``.
    """

    sample_id: str
    core_id: str                       # exact label in sdata.obs.core_id
    subject_id: str
    organ: str
    tma_slide: str
    condition: str
    tbd_status: str
    control_subtype: Optional[str]
    tma_sdata_path: Path               # shared across cores of the same TMA

    @property
    def core_obs_filter(self) -> Dict[str, str]:
        """Mask used by concat_samples to subset the TMA AnnData."""
        return {"core_id": self.core_id}


@dataclass
class Cohort:
    cohort_id: str
    description: str
    conditions: List[str]
    organs: Dict[str, Dict[str, Any]]
    donors: List[DonorSpec]
    integration: Dict[str, Any] = field(default_factory=dict)
    pseudobulk: Dict[str, Any] = field(default_factory=dict)
    ccc_multisample: Dict[str, Any] = field(default_factory=dict)
    spacia_meta: Dict[str, Any] = field(default_factory=dict)
    pseudotime_cohort: Dict[str, Any] = field(default_factory=dict)
    spatial_domains: Dict[str, Any] = field(default_factory=dict)

    # ─── selectors ──────────────────────────────────────────────────────
    def by_organ(self, organ: str) -> List[DonorSpec]:
        return [d for d in self.donors if d.organ == organ]

    def by_condition(self, condition: str) -> List[DonorSpec]:
        return [d for d in self.donors if d.condition == condition]

    def n_donors(self) -> int:
        return len(self.donors)

    def tma_sdata_for(self, organ: str) -> Path:
        if organ not in self.organs:
            raise KeyError(f"organ {organ!r} not in cohort.organs")
        return project_root() / self.organs[organ]["tma_sdata"]

    # ─── persistence ────────────────────────────────────────────────────
    def cohort_root(self) -> Path:
        return project_root() / "TBDs" / "cohort"

    def cohort_results(self) -> Path:
        return self.cohort_root() / "results"

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        return d

    def to_json(self, path: Path) -> None:
        """Dump as JSON for R-side scripts (muscat, slingshot)."""
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str))


# ────────────────────────────────────────────────────────────────────────
# Validation + parsing
# ────────────────────────────────────────────────────────────────────────
def _require(cfg: Dict[str, Any], key: str, where: str) -> Any:
    if key not in cfg:
        raise ValueError(f"{where}: missing required key '{key}'")
    return cfg[key]


def _parse_donor(raw: Dict[str, Any], cohort_conditions: List[str],
                 organs: Dict[str, Dict[str, Any]]) -> DonorSpec:
    where = f"donor '{raw.get('subject_id', '?')}'"
    sid = _require(raw, "subject_id", where)
    if not _SUBJECT_ID_RE.match(sid):
        raise ValueError(
            f"{where}: subject_id must match [A-Za-z0-9_]+ (got {sid!r})"
        )
    organ = _require(raw, "organ", where)
    if organ not in organs:
        raise ValueError(f"{where}: organ '{organ}' not declared in cohort.organs")
    tma_slide = _require(raw, "tma_slide", where)
    expected_tma = organs[organ].get("tma_slide")
    if expected_tma and tma_slide != expected_tma:
        raise ValueError(
            f"{where}: tma_slide '{tma_slide}' does not match organ default "
            f"'{expected_tma}'"
        )
    condition = _require(raw, "condition", where)
    if condition not in cohort_conditions:
        raise ValueError(
            f"{where}: condition '{condition}' not in cohort.conditions "
            f"{cohort_conditions}"
        )
    tbd = _require(raw, "tbd_status", where)
    if tbd not in _VALID_TBD:
        raise ValueError(
            f"{where}: tbd_status must be one of {_VALID_TBD} (got {tbd!r})"
        )
    return DonorSpec(
        subject_id=sid,
        organ=organ,
        tma_slide=tma_slide,
        tma_core_label=_require(raw, "tma_core_label", where),
        condition=condition,
        tbd_status=tbd,
        control_subtype=raw.get("control_subtype"),
        cores_expected=int(raw.get("cores_expected", 0)),
        notes=raw.get("notes", ""),
    )


def load_cohort(yaml_path: Optional[Path] = None) -> Cohort:
    path = Path(yaml_path) if yaml_path else project_root() / DEFAULT_COHORT_YAML
    if not path.exists():
        raise FileNotFoundError(f"cohort YAML not found: {path}")

    with open(path) as f:
        cfg = yaml.safe_load(f)

    where = f"cohort YAML {path.name}"
    cohort_id = _require(cfg, "cohort_id", where)
    conditions = _require(cfg, "conditions", where)
    if not isinstance(conditions, list) or len(conditions) < 2:
        raise ValueError(f"{where}: 'conditions' must be a list with ≥2 entries")
    organs = _require(cfg, "organs", where)
    if not isinstance(organs, dict) or not organs:
        raise ValueError(f"{where}: 'organs' must be a non-empty dict")
    for o, ocfg in organs.items():
        for k in ("tma_slide", "tma_sdata", "markers_yaml"):
            if k not in ocfg:
                raise ValueError(f"{where}: organs.{o} missing '{k}'")

    raw_donors = _require(cfg, "donors", where)
    if not raw_donors:
        raise ValueError(f"{where}: 'donors' list is empty")

    donors = [_parse_donor(d, conditions, organs) for d in raw_donors]
    seen: set[str] = set()
    for d in donors:
        if d.subject_id in seen:
            raise ValueError(f"{where}: duplicate subject_id {d.subject_id!r}")
        seen.add(d.subject_id)

    return Cohort(
        cohort_id=cohort_id,
        description=cfg.get("description", ""),
        conditions=conditions,
        organs=organs,
        donors=donors,
        integration=cfg.get("integration", {}),
        pseudobulk=cfg.get("pseudobulk", {}),
        ccc_multisample=cfg.get("ccc_multisample", {}),
        spacia_meta=cfg.get("spacia_meta", {}),
        pseudotime_cohort=cfg.get("pseudotime_cohort", {}),
        spatial_domains=cfg.get("spatial_domains", {}),
    )


# ────────────────────────────────────────────────────────────────────────
# Core expansion (TMA-aware)
# ────────────────────────────────────────────────────────────────────────
def expand_cores(
    cohort: Cohort,
    tma_obs: pd.DataFrame,
    tma_slide: str,
) -> List[SampleSpec]:
    """Expand donors into one SampleSpec per actual core on the TMA.

    Parameters
    ----------
    cohort
        Parsed cohort.
    tma_obs
        ``adata.obs`` of the per-TMA sdata.zarr (already dearrayed). Must
        contain columns ``core_id`` and ``tma_core_label``.
    tma_slide
        Slide identifier (e.g. 'LUNG_AF').

    Notes
    -----
    Matching is by ``tma_core_label``: every cell's label points to one
    donor, and every distinct ``core_id`` for that label becomes a sample.
    """
    required_cols = {"core_id", "tma_core_label"}
    missing = required_cols - set(tma_obs.columns)
    if missing:
        raise ValueError(
            f"tma_obs missing required columns: {missing}; "
            f"run TMA dearraying before calling expand_cores"
        )

    label_to_donor = {
        d.tma_core_label: d for d in cohort.donors if d.tma_slide == tma_slide
    }
    if not label_to_donor:
        raise ValueError(f"no donors registered for TMA slide {tma_slide!r}")

    # Unique (tma_core_label, core_id) pairs present on the slide
    pairs = tma_obs[list(required_cols)].drop_duplicates().reset_index(drop=True)

    samples: List[SampleSpec] = []
    counts: Dict[str, int] = {d.subject_id: 0 for d in label_to_donor.values()}
    organs_lookup = {d.tma_slide: d.organ for d in cohort.donors}

    for _, row in pairs.iterrows():
        label = row["tma_core_label"]
        core_id = row["core_id"]
        donor = label_to_donor.get(label)
        if donor is None:
            continue                  # label not in cohort manifest — skip
        counts[donor.subject_id] += 1
        idx = counts[donor.subject_id]
        sample_id = f"{donor.subject_id}_core{idx:02d}"
        samples.append(SampleSpec(
            sample_id=sample_id,
            core_id=str(core_id),
            subject_id=donor.subject_id,
            organ=donor.organ,
            tma_slide=donor.tma_slide,
            condition=donor.condition,
            tbd_status=donor.tbd_status,
            control_subtype=donor.control_subtype,
            tma_sdata_path=cohort.tma_sdata_for(donor.organ),
        ))

    return samples


def cohort_summary(cohort: Cohort) -> pd.DataFrame:
    """Donor-level summary table for QC/reporting."""
    rows = []
    for d in cohort.donors:
        rows.append({
            "subject_id": d.subject_id,
            "organ": d.organ,
            "tma_slide": d.tma_slide,
            "condition": d.condition,
            "tbd_status": d.tbd_status,
            "control_subtype": d.control_subtype or "",
            "cores_expected": d.cores_expected,
        })
    return pd.DataFrame(rows)
