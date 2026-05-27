#!/usr/bin/env bash
# Compute PGS from the bundled score files using PLINK2.
#
# Edit the variables below for your environment, then run:
#   bash 01_compute_pgs.sh
#
# Score file columns: 1=SNP, 2=A1, 3=A2, 4..N=PGS columns.
# Output for each score file: <out_prefix>.sscore with one column per PGS column.

set -euo pipefail

# -------- EDIT THESE --------
PLINK2=plink2                                                # path to PLINK2 binary
TARGET_PFILE=/path/to/your/target_data                       # PLINK2 .pgen/.pvar/.psam prefix
SCORES_DIR=../scores                                         # relative to this script
OUT_DIR=../pgs_out
# ----------------------------

mkdir -p "$OUT_DIR"

compute() {
    local score_file=$1
    local out_prefix=$2

    # Number of score columns = total columns - 3 (SNP, A1, A2)
    local ncol
    ncol=$(zcat "$score_file" | head -1 | awk '{print NF}')
    local last=$ncol

    echo ">>> Scoring with $(basename "$score_file"): $((ncol - 3)) PGS columns"
    "$PLINK2" \
        --pfile "$TARGET_PFILE" \
        --score "$score_file" 1 2 header-read cols=+scoresums \
        --score-col-nums 4-"$last" \
        --out "$out_prefix"
}

compute "$SCORES_DIR/pseudovalidated.score.gz" "$OUT_DIR/pseudovalidated"

# Uncomment to also score the full hyperparameter grid (much larger):
# compute "$SCORES_DIR/full_grid.score.gz" "$OUT_DIR/full_grid"

echo "Done. Output in: $OUT_DIR"
echo "Inspect e.g.: head $OUT_DIR/pseudovalidated.sscore"
echo
echo "Report PLINK2's variant count line (e.g. '--score: X variants processed.')"
echo "to confirm overlap with the target sample is high."
