#!/bin/bash
# Round 7b Section 4 driver — run LDpred2 pipeline for each of the 6 candidate
# panels in parallel waves.
#
# Prerequisites:
#   - setup_r7bS4_ldpred2.sh already run (6 private ldpred2_ldref dirs + 6 configs)
#   - Container active with genopred env
#
# Usage:
#   bash run_r7bS4_ldpred2.sh <wave>
# where <wave> ∈ {1, 2, all}. Wave 1: MELD_lambda0, EUR, EAS. Wave 2: AFR, CSA, AMR.

set -uo pipefail

if [ $# -lt 1 ]; then echo "usage: $0 <wave>"; exit 1; fi
WAVE=$1
BASE=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline

case "$WAVE" in
  1)   LABELS="MELD_lambda0 EUR EAS" ;;
  2)   LABELS="AFR CSA AMR" ;;
  all) LABELS="MELD_lambda0 EUR EAS AFR CSA AMR" ;;
  *) echo "wave must be 1, 2, or all"; exit 1 ;;
esac

cd "$BASE"

# Ensure snakemake and env are active — genopred env should already be sourced.
if ! command -v snakemake >/dev/null 2>&1; then
  # Fall back to direct path
  export PATH=/home/claude/micromamba/envs/genopred/bin:$PATH
fi

mkdir -p /home/claude/r7bS4_ldpred2_logs

echo "=== Wave $WAVE — labels: $LABELS ==="
declare -A PIDS
for L in $LABELS; do
  CFG=misc/opensnp/config_meld_r7bS4_ldpred2_${L}.yaml
  LOG=/home/claude/r7bS4_ldpred2_logs/${L}.log
  echo "launching $L → $LOG"
  nohup snakemake --use-conda --conda-frontend conda -j 5 \
    --configfile=$CFG \
    prep_pgs_ldpred2 \
    > "$LOG" 2>&1 &
  PIDS[$L]=$!
done

# Wait for all
FAIL=0
for L in $LABELS; do
  PID=${PIDS[$L]}
  wait "$PID" && echo "$L: OK" || { echo "$L: FAILED (exit $?)"; FAIL=1; }
done

exit $FAIL
