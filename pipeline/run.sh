#!/usr/bin/env bash
# =============================================================================
# run.sh — Xenium pipeline entry point
# Usage:
#   ./run.sh                          # full pipeline
#   ./run.sh --until ingest           # run only up to ingest
#   ./run.sh --cores 8                # override core count
#   ./run.sh --forcerun qc            # re-run a specific step
#   ./run.sh --configfile config/config_lung.yaml --until ingest
# =============================================================================

set -euo pipefail

# ─── Defaults ─────────────────────────────────────────────────────────────────
CONFIGFILE="config/config_lung.yaml"
CORES=8
EXTRA_ARGS=()

# ─── Parse arguments ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --configfile) CONFIGFILE="$2"; shift 2 ;;
        --cores)      CORES="$2";      shift 2 ;;
        *)            EXTRA_ARGS+=("$1"); shift ;;
    esac
done

# ─── Setup ────────────────────────────────────────────────────────────────────
RESULTS_DIR=$(python3 -c "
import yaml
with open('$CONFIGFILE') as f:
    cfg = yaml.safe_load(f)
print(cfg['project']['output_dir'])
")
LOGS_DIR="${RESULTS_DIR}/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MAIN_LOG="${LOGS_DIR}/pipeline_${TIMESTAMP}.log"
SUMMARY="${RESULTS_DIR}/pipeline_summary.txt"
STATS_JSON="${LOGS_DIR}/pipeline_stats_${TIMESTAMP}.json"

mkdir -p "$LOGS_DIR"

echo "============================================================"
echo " Xenium Pipeline"
echo " Config:    $CONFIGFILE"
echo " Cores:     $CORES"
echo " Log:       $MAIN_LOG"
echo " Started:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"

# ─── Run Snakemake ────────────────────────────────────────────────────────────
START_TIME=$(date +%s)

snakemake \
    --configfile "$CONFIGFILE" \
    --cores "$CORES" \
    --stats "$STATS_JSON" \
    --printshellcmds \
    --show-failed-logs \
    --rerun-incomplete \
    "${EXTRA_ARGS[@]}" \
    2>&1 | tee "$MAIN_LOG"

EXIT_CODE=${PIPESTATUS[0]}
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))

# ─── Generate summary ─────────────────────────────────────────────────────────
python3 scripts/utils/summarize_run.py \
    --logs-dir      "$LOGS_DIR" \
    --stats-json    "$STATS_JSON" \
    --main-log      "$MAIN_LOG" \
    --output        "$SUMMARY" \
    --exit-code     "$EXIT_CODE" \
    --elapsed       "$ELAPSED" \
    --configfile    "$CONFIGFILE"

echo ""
echo "============================================================"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo " PIPELINE COMPLETED SUCCESSFULLY"
else
    echo " PIPELINE FAILED (exit code: $EXIT_CODE)"
fi
echo " Elapsed:   $(printf '%02dh %02dm %02ds' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))"
echo " Summary:   $SUMMARY"
echo " Full log:  $MAIN_LOG"
echo "============================================================"

exit $EXIT_CODE
