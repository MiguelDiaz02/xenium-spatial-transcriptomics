#!/usr/bin/env bash
# Phase F — Run all 30 Spacia jobs sequentially (avoids RAM saturation with 268k cells).
# Each job takes ~5-15 min with -m 20000,10000,50,2 -nc 3000. Total ~4-6 hours.
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
INTERACTIONS="$DATA_DIR/F3_selected_interactions.csv"

if [[ ! -f "$SPACIA" ]]; then
    echo "ERROR: spacia.py not found at $SPACIA" >&2
    echo "Set SPACIA_PATH to the location of spacia.py or place it under pipeline/external/Spacia/" >&2
    exit 1
fi

mkdir -p "$OUTBASE"

MCMC="20000,10000,50,2"
N_CELLS=3000
DIST=30

SUCCESS=0
FAIL=0

# Read interactions CSV (skip header)
tail -n +2 "$INTERACTIONS" | while IFS=',' read -r source target ligand receptor rest; do
    # Strip quotes if present
    source=$(echo "$source" | tr -d '"')
    target=$(echo "$target" | tr -d '"')
    ligand=$(echo "$ligand" | tr -d '"')
    receptor=$(echo "$receptor" | tr -d '"')

    # Use first gene for complex ligands/receptors
    sf=$(echo "$ligand" | cut -d'|' -f1)
    rf=$(echo "$receptor" | cut -d'|' -f1)

    JOB_NAME="${source}__${target}__${ligand}__${receptor}"
    JOB_OUT="$OUTBASE/$JOB_NAME"

    # Skip if already done (B_and_FDR.csv is the success marker)
    if [ -f "$JOB_OUT/B_and_FDR.csv" ]; then
        echo "[CACHED] $JOB_NAME"
        continue
    fi

    mkdir -p "$JOB_OUT"
    echo "[RUN] $source -> $target | $sf -> $rf"

    conda run -n spacia python "$SPACIA" "$COUNTS" "$META" \
        -sc "$source" -rc "$target" \
        -sf "$sf" -rf "$rf" \
        -d $DIST -m $MCMC -nc $N_CELLS \
        -o "$JOB_OUT" \
        > "$JOB_OUT/stdout.log" 2>&1

    if [ -f "$JOB_OUT/B_and_FDR.csv" ]; then
        echo "[OK] $JOB_NAME"
        SUCCESS=$((SUCCESS+1))
    else
        echo "[FAIL] $JOB_NAME"
        FAIL=$((FAIL+1))
    fi
done

echo "Done: $SUCCESS succeeded, $FAIL failed"
