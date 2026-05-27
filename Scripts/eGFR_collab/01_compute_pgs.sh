#!/usr/bin/env bash
# Compute PGS from the bundled score files using PLINK2.
#
# Two modes:
#   * Single-prefix target data: set TARGET_PFILE to a single .pgen/.pvar/.psam prefix.
#   * Per-chromosome target data: set TARGET_PFILE_CHR to a prefix with {CHR}
#     placeholder, e.g. /data/chr{CHR}/genotypes . The script then scores each
#     chromosome separately and sums the per-chromosome PGS sums to produce
#     a single .sscore per score file.
#
# Score file columns: 1=SNP, 2=A1, 3=A2, 4..N=PGS columns.
# Output: <OUT_DIR>/<basename>.sscore with one _AVG and one _SUM column per PGS.

set -euo pipefail

# -------- EDIT THESE --------
PLINK2=plink2                                                # path to PLINK2 binary
TARGET_PFILE=                                                # e.g. /path/to/all_chr_data  (leave blank if per-chr)
TARGET_PFILE_CHR=                                            # e.g. /path/to/chr{CHR}/genotypes (use {CHR})
CHROMOSOMES=$(seq 1 22)                                      # space-separated chromosome list
SCORES_DIR=../scores                                         # relative to this script
OUT_DIR=../pgs_out
SCORE_FILES=("$SCORES_DIR/pseudovalidated.score.gz")         # add full_grid.score.gz here if needed
# ----------------------------

mkdir -p "$OUT_DIR"

score_single() {
    local score_file=$1
    local out_prefix=$2
    local pfile=$3

    local ncol
    ncol=$(zcat "$score_file" | head -1 | awk '{print NF}')

    "$PLINK2" \
        --pfile "$pfile" \
        --score "$score_file" 1 2 header-read cols=+scoresums \
        --score-col-nums 4-"$ncol" \
        --out "$out_prefix"
}

score_perchr_and_sum() {
    # Score each chromosome separately, then sum *_SUM columns and recompute *_AVG
    # over the union of allele counts. Uses awk for the final aggregation.
    local score_file=$1
    local out_prefix=$2
    local tmpdir
    tmpdir=$(mktemp -d)

    for CHR in $CHROMOSOMES; do
        local pfile="${TARGET_PFILE_CHR//\{CHR\}/$CHR}"
        [[ -e "${pfile}.pgen" ]] || { echo "Skip chr$CHR — ${pfile}.pgen not found"; continue; }
        score_single "$score_file" "$tmpdir/chr$CHR" "$pfile" >/dev/null
    done

    # Combine using R: stack per-chr sscore files, sum summable cols by sample ID.
    Rscript --vanilla - "$tmpdir" "$out_prefix.sscore" <<'EOF'
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
tmpdir <- args[1]; out <- args[2]
files <- list.files(tmpdir, pattern = "\\.sscore$", full.names = TRUE)
if (length(files) == 0) stop("No per-chromosome .sscore files produced")
stacked <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
nm <- names(stacked)
id_col   <- intersect(c("#IID", "IID"), nm)[1]
sum_cols <- unique(c(intersect(c("ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM"), nm),
                     grep("_SUM$", nm, value = TRUE)))
combined <- stacked[, lapply(.SD, sum, na.rm = TRUE), by = id_col, .SDcols = sum_cols]
# Recompute _AVG = _SUM / ALLELE_CT (skip the bookkeeping NAMED_ALLELE_DOSAGE_SUM)
score_sum_cols <- setdiff(grep("_SUM$", names(combined), value = TRUE),
                          "NAMED_ALLELE_DOSAGE_SUM")
for (s in score_sum_cols) {
  combined[[sub("_SUM$", "_AVG", s)]] <- combined[[s]] / combined$ALLELE_CT
}
fwrite(combined, out, sep = "\t")
cat(sprintf("Wrote %s (%d samples)\n", out, nrow(combined)))
EOF
    rm -rf "$tmpdir"
}

run() {
    local score_file=$1
    local out_prefix="$OUT_DIR/$(basename "$score_file" .score.gz)"
    echo ">>> Scoring with $(basename "$score_file")"
    if [[ -n "$TARGET_PFILE" ]]; then
        score_single "$score_file" "$out_prefix" "$TARGET_PFILE"
    elif [[ -n "$TARGET_PFILE_CHR" ]]; then
        score_perchr_and_sum "$score_file" "$out_prefix"
    else
        echo "ERROR: set either TARGET_PFILE or TARGET_PFILE_CHR" >&2
        exit 1
    fi
}

for sf in "${SCORE_FILES[@]}"; do
    [[ -e "$sf" ]] && run "$sf" || echo "Missing: $sf — skipping"
done

echo
echo "Done. Output in: $OUT_DIR"
echo
echo "Sanity check: PLINK2 reports '--score: N variants processed.' per run."
echo "Sum across chromosomes should be close to the total in the score file"
echo "(~1.2M for this bundle)."
