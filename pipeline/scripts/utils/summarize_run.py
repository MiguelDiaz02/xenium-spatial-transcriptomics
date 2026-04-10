"""
Post-run summary generator.
Reads per-step logs + Snakemake stats JSON → pipeline_summary.txt

Output is plain text designed to be:
  - Human-readable at a glance
  - Grep-friendly for errors
  - Compact (no redundant info)
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import timedelta
from pathlib import Path


STEP_ORDER = [
    "01_ingest",
    "02_segmentation",
    "02b_compare_segmentation",
    "03_qc",
    "04_preprocess",
    "05_reduction",
    "06_annotation",
    "07_denoising",
    "08_spatial",
    "09_downstream",
]

# Patterns to extract key metrics from each step log
METRIC_PATTERNS = {
    "01_ingest": [
        (r"SpatialData loaded", "SpatialData created"),
        (r"Write complete", "Zarr store written"),
        (r"H&E image registered", "H&E registered"),
        (r"H&E registration failed.*?;", "WARNING: H&E registration failed"),
    ],
    "02_segmentation": [
        (r"Cellpose detected ([\d,]+) cells", "Cells detected (Cellpose): {}"),
        (r"Stored ([\d,]+) Cellpose polygons", "Polygons stored: {}"),
        (r"Stored ([\d,]+) Baysor polygons",  "Polygons stored (Baysor): {}"),
        (r"method=xoa.*no changes", "Method: XOA (no re-segmentation)"),
    ],
    "03_qc": [
        (r"QC filter: ([\d,]+) → ([\d,]+) cells \((.+?) removed", "Cells: {} → {} ({} removed)"),
        (r"Found (\d+) negative control genes", "Neg ctrl genes: {}"),
    ],
    "04_preprocess": [
        (r"Saving raw counts", "Raw counts saved to layers['counts']"),
        (r"Normalizing total counts \(target_sum=(\S+)\)", "Normalized (target_sum={})"),
        (r"After HVG selection: (\d+) genes", "HVG selected: {} genes"),
        (r"Skipping HVG", "HVG selection: skipped (small panel)"),
    ],
    "05_reduction": [
        (r"PCA: (\d+) components on ([\d,]+) cells", "PCA: {} components on {} cells"),
        (r"Leiden: (\d+) clusters detected", "Clusters (Leiden): {}"),
    ],
    "06_annotation": [
        (r"CellTypist: (\d+) cell types annotated", "Cell types (CellTypist): {}"),
        (r"Manual annotation: (\d+) cell types", "Cell types (manual): {}"),
    ],
    "07_denoising": [
        (r"Training ResolVI", "ResolVI training started"),
        (r"Denoised layer added: shape (.+)", "Denoised layer shape: {}"),
        (r"Raw counts layer preserved: ✓", "Raw counts: preserved ✓"),
        (r"Skipped.*denoising_resolvi", "ResolVI: skipped (disabled in config)"),
    ],
    "08_spatial": [
        (r"Top spatially variable genes", "Moran's I computed"),
        (r"Neighborhood enrichment", "Neighborhood enrichment computed"),
        (r"Co-occurrence analysis", "Co-occurrence computed"),
        (r"Figures saved to (.+)", "Figures: {}"),
    ],
    "09_downstream": [
        (r"Scrublet: ([\d,]+) predicted doublets \((.+?)\)", "Doublets: {} ({})"),
        (r"(\d+) significant DEGs", "Significant DEGs: {}"),
        (r"Pseudobulk DE: skipped", "Pseudobulk: skipped (single sample)"),
    ],
}


def parse_step_log(log_path: Path) -> dict:
    """Extract status, metrics, errors, and warnings from a step log."""
    result = {
        "status": "not_run",
        "metrics": [],
        "errors": [],
        "warnings": [],
    }

    if not log_path.exists():
        return result

    text = log_path.read_text(errors="replace")
    lines = text.splitlines()

    # Determine status
    step_name = log_path.stem  # e.g. "01_ingest"
    done_keywords = ["DONE", "complete", "Complete", "successfully"]
    fail_keywords  = ["Error", "ERROR", "Traceback", "Exception", "FAILED"]

    for line in reversed(lines):
        if any(k in line for k in done_keywords):
            result["status"] = "success"
            break
        if any(k in line for k in fail_keywords):
            result["status"] = "failed"
            break
    else:
        result["status"] = "incomplete"

    # Extract errors with context
    for i, line in enumerate(lines):
        if any(k in line for k in ["ERROR", "Traceback", "Exception"]):
            # Include up to 3 lines of context
            context = lines[max(0, i-1):min(len(lines), i+4)]
            result["errors"].append("\n    ".join(c.strip() for c in context if c.strip()))
        elif "WARNING" in line or "warning" in line:
            result["warnings"].append(line.strip())

    # Extract metrics
    patterns = METRIC_PATTERNS.get(step_name, [])
    for pattern, template in patterns:
        match = re.search(pattern, text)
        if match:
            groups = match.groups()
            if groups:
                try:
                    metric = template.format(*groups)
                except (IndexError, KeyError):
                    metric = template
            else:
                metric = template
            result["metrics"].append(metric)

    return result


def parse_stats_json(stats_path: Path) -> dict[str, float]:
    """Extract per-rule runtime from Snakemake stats JSON."""
    if not stats_path or not stats_path.exists():
        return {}
    try:
        data = json.loads(stats_path.read_text())
        rules = data.get("rules", {})
        return {rule: info.get("mean-runtime", 0) for rule, info in rules.items()}
    except Exception:
        return {}


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f}s"
    td = timedelta(seconds=int(seconds))
    h, rem = divmod(td.seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h {m:02d}m {s:02d}s"
    return f"{m}m {s:02d}s"


def scan_main_log_for_errors(main_log: Path) -> list[str]:
    """Pull error blocks from the main Snakemake log."""
    if not main_log.exists():
        return []
    errors = []
    lines = main_log.read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if "Error in rule" in line or "RuleException" in line or "WorkflowError" in line:
            block = lines[i:min(len(lines), i+10)]
            errors.append("\n  ".join(b.strip() for b in block if b.strip()))
    return errors


def build_summary(
    logs_dir: Path,
    stats_json: Path,
    main_log: Path,
    exit_code: int,
    elapsed: int,
    configfile: str,
) -> str:
    runtimes = parse_stats_json(stats_json)
    main_errors = scan_main_log_for_errors(main_log)

    lines = []
    sep = "=" * 62

    lines += [
        sep,
        "  XENIUM PIPELINE — RUN SUMMARY",
        sep,
        f"  Config:   {configfile}",
        f"  Log dir:  {logs_dir}",
        f"  Elapsed:  {format_elapsed(elapsed)}",
        f"  Status:   {'SUCCESS ✓' if exit_code == 0 else f'FAILED ✗ (exit {exit_code})'}",
        sep,
        "",
    ]

    all_errors   = []
    all_warnings = []

    lines.append("  STEPS")
    lines.append("  " + "-" * 58)

    for step in STEP_ORDER:
        log_path = logs_dir / f"{step}.log"
        info = parse_step_log(log_path)

        status_icon = {"success": "✓", "failed": "✗", "not_run": "—", "incomplete": "?"}.get(
            info["status"], "?"
        )

        # Map step name to Snakemake rule name for runtime lookup
        rule_map = {
            "01_ingest": "ingest",
            "02_segmentation": "segmentation",
            "02b_compare_segmentation": "compare_segmentation",
            "03_qc": "qc",
            "04_preprocess": "preprocess",
            "05_reduction": "reduction",
            "06_annotation": "annotation",
            "07_denoising": "denoising",
            "08_spatial": "spatial",
            "09_downstream": "downstream",
        }
        rule = rule_map.get(step, step)
        rt = runtimes.get(rule, 0)
        rt_str = f"  {format_elapsed(rt):>8}" if rt else "          "

        lines.append(f"  {status_icon} {step:<32}{rt_str}")

        for metric in info["metrics"]:
            lines.append(f"      → {metric}")

        if info["warnings"]:
            for w in info["warnings"][:3]:   # cap at 3 warnings per step
                lines.append(f"      ⚠ {w[:80]}")
            all_warnings.extend(info["warnings"])

        if info["errors"]:
            for e in info["errors"]:
                lines.append(f"      ✗ {e[:120]}")
            all_errors.extend(info["errors"])

    lines.append("")

    # Errors section
    if all_errors or main_errors:
        lines += [sep, "  ERRORS", "  " + "-" * 58]
        for err in (all_errors + main_errors):
            lines.append(f"  {err[:200]}")
            lines.append("")
    else:
        lines += [sep, "  ERRORS:   None", ""]

    # Warnings section
    if all_warnings:
        lines += [sep, "  WARNINGS", "  " + "-" * 58]
        for w in all_warnings:
            lines.append(f"  {w[:120]}")
    else:
        lines += [sep, "  WARNINGS: None", ""]

    lines += ["", sep]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--logs-dir",   required=True)
    parser.add_argument("--stats-json", required=False, default="")
    parser.add_argument("--main-log",   required=True)
    parser.add_argument("--output",     required=True)
    parser.add_argument("--exit-code",  type=int, default=0)
    parser.add_argument("--elapsed",    type=int, default=0)
    parser.add_argument("--configfile", default="")
    args = parser.parse_args()

    summary = build_summary(
        logs_dir   = Path(args.logs_dir),
        stats_json = Path(args.stats_json) if args.stats_json else Path(""),
        main_log   = Path(args.main_log),
        exit_code  = args.exit_code,
        elapsed    = args.elapsed,
        configfile = args.configfile,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
