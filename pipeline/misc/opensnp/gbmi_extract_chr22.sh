#!/usr/bin/env bash
# Extract chr22 from every downloaded GBMI .gz sumstat into a per-file
# .chr22.tsv.gz. Uses tabix if the .tbi is present (fast, seek-based) and
# falls back to zcat+awk otherwise (slower, full scan) so we can start
# processing before the whole download completes.
#
# Output: alongside the source .gz, with .chr22.tsv.gz suffix.
# The extracted files keep the original 17-column header prefixed with `#CHR`.

set -eu -o pipefail

DEST=/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI
OUT=$DEST/chr22
LOG=$DEST/gbmi_r8_extract.log
TABIX=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/tabix
BGZIP=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/bgzip

mkdir -p "$OUT"

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) start extract ===" | tee -a "$LOG"

usage_tabix=0
usage_awk=0

for src in "$DEST"/*_nbbkgt1.txt.gz; do
  fn=$(basename "$src")
  out=$OUT/${fn%.txt.gz}.chr22.tsv.gz

  if [[ -f "$out" && -s "$out" ]]; then
    echo "skip (already exists): $out" | tee -a "$LOG"
    continue
  fi

  if [[ -f "${src}.tbi" ]]; then
    # tabix seek — fast
    "$TABIX" -h "$src" 22 | "$BGZIP" -c > "$out"
    usage_tabix=$((usage_tabix+1))
    method=tabix
  else
    # awk scan — works without .tbi
    { zcat "$src" | head -1; zcat "$src" | awk 'BEGIN{FS="\t"; OFS="\t"} $1=="22"'; } \
      | "$BGZIP" -c > "$out"
    usage_awk=$((usage_awk+1))
    method=awk
  fi

  rows=$(zcat "$out" | wc -l)
  size=$(du -h "$out" | cut -f1)
  echo "[$method] $fn -> $rows rows, $size" | tee -a "$LOG"
done

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) done: tabix=$usage_tabix awk=$usage_awk ===" | tee -a "$LOG"
