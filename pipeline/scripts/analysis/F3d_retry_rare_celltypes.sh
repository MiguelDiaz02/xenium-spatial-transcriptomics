#!/usr/bin/env bash
# Retry Spacia jobs that failed due to rare sender cell types (<500 cells).
# Uses -nc 100 (small subsample) instead of -nc 3000.
#
# Path resolution:
#   XENIUM_PROJECT_ROOT  — project root (auto-detected from this script's path).
#   XENIUM_DATASET       — dataset name under PROJECT_ROOT (default: human_lung_cancer).
#   SPACIA_PATH          — path to spacia.py (default: $PROJECT_ROOT/pipeline/external/Spacia/spacia.py).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${XENIUM_PROJECT_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}"
DATASET="${XENIUM_DATASET:-human_lung_cancer}"
SPACIA="${SPACIA_PATH:-$PROJECT_ROOT/pipeline/external/Spacia/spacia.py}"

DATA_DIR="$PROJECT_ROOT/$DATASET/results/02_biology/phase_f_spacia"
COUNTS="$DATA_DIR/spacia_counts.txt"
META="$DATA_DIR/spacia_metadata.txt"
OUTBASE="$DATA_DIR/jobs"

if [[ ! -f "$SPACIA" ]]; then
    echo "ERROR: spacia.py not found at $SPACIA" >&2
    echo "Set SPACIA_PATH to the location of spacia.py or place it under pipeline/external/Spacia/" >&2
    exit 1
fi

MCMC="20000,10000,50,2"
N_CELLS=100   # rare cell types: subsample conservatively
DIST=30

SUCCESS=0; FAIL=0

run_job() {
    local source=$1 target=$2 sf=$3 rf=$4
    local JOB_NAME="${source}__${target}__${sf}__${rf}"
    local JOB_OUT="$OUTBASE/$JOB_NAME"

    if [ -f "$JOB_OUT/B_and_FDR.csv" ]; then
        echo "[CACHED] $JOB_NAME"; return
    fi

    rm -rf "$JOB_OUT" && mkdir -p "$JOB_OUT"
    echo "[RETRY -nc $N_CELLS] $source -> $target | $sf -> $rf"

    conda run -n spacia python "$SPACIA" "$COUNTS" "$META" \
        -sc "$source" -rc "$target" \
        -sf "$sf" -rf "$rf" \
        -d $DIST -m $MCMC -nc $N_CELLS \
        -o "$JOB_OUT" \
        > "$JOB_OUT/stdout.log" 2>&1

    if [ -f "$JOB_OUT/B_and_FDR.csv" ]; then
        echo "[OK] $JOB_NAME"; SUCCESS=$((SUCCESS+1))
    else
        echo "[FAIL] $JOB_NAME (rare cell type — insufficient bags)"; FAIL=$((FAIL+1))
    fi
}

# Known failures: rare sender types
run_job "Plasma" "Tumor_resting" "ADAM17" "MUC1"
run_job "pDC" "Tumor_resting" "ADAM17" "MUC1"
run_job "Plasma" "Epithelial_general" "ADAM17" "MUC1"
run_job "pDC" "Epithelial_general" "ADAM17" "MUC1"
N_CELLS=50  # DC_mature has only ~98 cells
run_job "DC_mature" "Tumor_resting" "ADAM17" "MUC1"
N_CELLS=100
run_job "Plasmablast" "Tumor_resting" "ADAM17" "MUC1"

echo "Retry done: $SUCCESS OK, $FAIL failed"
