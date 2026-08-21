#!/usr/bin/env Rscript
# Round 7b Section 4 M5 scorer — pairwise weight-vector correlations across LDpred2
# panels for the same Yengo-2022 all-ancestry GWAS. Intersects SNP sets, computes
# per-model pairwise Pearson correlations, block-bootstrap CIs.
#
# Input: six per-panel LDpred2 score files, one per candidate label. Score files
# are gzipped tables at
#   /users/.../GenoPred/meld_test_r7bS4_ldpred2_<label>/reference/pgs_score_files/ldpred2/yengo_all/ref-yengo_all.score.gz
#
# Output (under pipeline/misc/opensnp/):
#   meld_r7b_S4_m5_intersect_summary.csv    per-panel drop counts + intersection size
#   meld_r7b_S4_m5_pairwise.csv             pairwise per-model correlations with bootstrap CIs
#   meld_r7b_S4_m5_pairwise_matrix.csv      symmetric matrix for the primary model
#
# Deliverable table on stdout for the report.

suppressPackageStartupMessages({
  library(data.table)
})

BASE_OUT <- '/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred'
LABELS   <- c('MELD_lambda0','EUR','EAS','AFR','CSA','AMR')
OUT_DIR  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
B_BOOT   <- 2000L
set.seed(2026L)

sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# ---------------------------------------------------------------------------
# Load per-panel score files
score_paths <- setNames(
  file.path(BASE_OUT,
            sprintf('meld_test_r7bS4_ldpred2_%s', LABELS),
            'reference/pgs_score_files/ldpred2/yengo_all/ref-yengo_all.score.gz'),
  LABELS)

panels <- list()
for (L in LABELS) {
  fp <- score_paths[[L]]
  if (!file.exists(fp)) { cat(sprintf('MISSING: %s (%s) — skipping\n', L, fp)); next }
  d <- fread(fp)
  cat(sprintf('%s: %d rows, cols: %s\n', L, nrow(d), paste(names(d), collapse = ', ')))
  panels[[L]] <- d
}
if (length(panels) < 2L) stop('need >= 2 panel score files')

# Detect the SNP identifier and effect columns present.
# ldpred2.R score files typically have SNP, A1, A2, and a set of score columns
# like SCORE_ldpred2_auto, SCORE_ldpred2_inf, SCORE_ldpred2_grid_*
snp_col_candidates <- c('SNP','ID','rsid','variant_id')
snp_col <- intersect(snp_col_candidates, names(panels[[1L]]))[1L]
if (is.na(snp_col)) stop('cannot find SNP column in ', LABELS[1L])
cat(sprintf('SNP identifier column: %s\n', snp_col))

# Find score / effect columns (starting with 'SCORE')
score_cols <- grep('^SCORE', names(panels[[1L]]), value = TRUE)
if (length(score_cols) == 0L) {
  # Fallback: any numeric column not obviously metadata
  score_cols <- setdiff(names(panels[[1L]]),
                        c(snp_col,'CHR','BP','A1','A2','SNP','ID','POS','pos','chr'))
  score_cols <- score_cols[sapply(score_cols, function(c) is.numeric(panels[[1L]][[c]]))]
}
cat(sprintf('score columns: %s\n', paste(score_cols, collapse = ', ')))

# ---------------------------------------------------------------------------
# Intersect SNPs
snp_sets <- lapply(panels, function(d) d[[snp_col]])
intersection <- Reduce(intersect, snp_sets)
cat(sprintf('intersection size: %d SNPs\n', length(intersection)))

intersect_summary <- rbindlist(lapply(names(panels), function(L) {
  n0 <- length(snp_sets[[L]])
  ni <- length(intersection)
  data.table(panel = L, n_panel = n0, n_intersect = ni,
             n_dropped = n0 - ni, frac_kept = ni / n0)
}))
fwrite(intersect_summary, file.path(OUT_DIR, 'meld_r7b_S4_m5_intersect_summary.csv'))

# Aligned per-panel effect vectors on intersection, per score column.
# Also align alleles: if a panel's A1 for a given SNP disagrees with the reference
# panel (first label), flip the effect sign. The score file A1 should be the
# effect-coded allele.
ref_panel <- panels[[LABELS[1L]]]
ref_a1 <- setNames(ref_panel$A1[match(intersection, ref_panel[[snp_col]])],
                   intersection)

align_and_extract <- function(d, sc) {
  ord <- match(intersection, d[[snp_col]])
  eff <- d[[sc]][ord]
  a1  <- d$A1[ord]
  sign <- ifelse(a1 == ref_a1, 1, ifelse(!is.na(a1) & a1 != ref_a1, -1, NA))
  eff * sign
}

# ---------------------------------------------------------------------------
# Pairwise correlations per score column, with block-bootstrap CIs.
# For block-bootstrap: group SNPs by LDetect block using the ldpred2 map.
# Fallback if map isn't available: use a large-block chunking heuristic (chr).
group_id_by_snp <- NULL
map_path <- file.path(OUT_DIR, 'ldpred2_ref_meld', LABELS[1L], 'map.rds')
if (file.exists(map_path)) {
  m <- readRDS(map_path)
  group_id_by_snp <- setNames(m$group_id, m$rsid)
  # Global unique block id = chr * 10000 + group_id (block_id was per-chr)
  chr_by_snp <- setNames(m$chr, m$rsid)
  gb <- 10000L * chr_by_snp[intersection] + group_id_by_snp[intersection]
} else {
  cat('WARNING: no map.rds found, using chr for block bootstrap grouping\n')
  gb <- rep(NA_integer_, length(intersection))
}
block_ids <- unique(gb[!is.na(gb)])
cat(sprintf('block-bootstrap group count: %d\n', length(block_ids)))

boot_paired_r <- function(a, b, gb, B = B_BOOT) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5L) return(c(mean = NA_real_, sd = NA_real_,
                             lo95 = NA_real_, hi95 = NA_real_))
  a <- a[ok]; b <- b[ok]; gb <- gb[ok]
  r0 <- suppressWarnings(cor(a, b))
  bids <- unique(gb[!is.na(gb)])
  # Fall back to SNP-level bootstrap if no block structure
  if (length(bids) < 5L) {
    boots <- replicate(B, {
      i <- sample.int(length(a), length(a), replace = TRUE)
      suppressWarnings(cor(a[i], b[i]))
    })
  } else {
    # index into rows per block
    row_by_block <- split(seq_along(a), gb)
    row_by_block <- row_by_block[names(row_by_block) != 'NA']
    bnames <- names(row_by_block)
    boots <- replicate(B, {
      pick <- sample(bnames, length(bnames), replace = TRUE)
      idx  <- unlist(row_by_block[pick])
      suppressWarnings(cor(a[idx], b[idx]))
    })
  }
  c(mean = r0, sd = sd(boots, na.rm = TRUE),
    lo95 = as.numeric(quantile(boots, 0.025, na.rm = TRUE)),
    hi95 = as.numeric(quantile(boots, 0.975, na.rm = TRUE)))
}

pair_rows <- list()
for (sc in score_cols) {
  # Extract aligned vectors for every panel
  vecs <- lapply(panels, align_and_extract, sc = sc)
  panel_names <- names(vecs)
  for (i in 1:(length(panel_names) - 1L)) {
    for (j in (i + 1L):length(panel_names)) {
      st <- boot_paired_r(vecs[[panel_names[i]]], vecs[[panel_names[j]]], gb)
      pair_rows[[length(pair_rows) + 1L]] <- data.table(
        model = sc, panel_a = panel_names[i], panel_b = panel_names[j],
        r = st[['mean']], sd = st[['sd']],
        lo95 = st[['lo95']], hi95 = st[['hi95']])
    }
  }
}
pair_dt <- rbindlist(pair_rows)
fwrite(pair_dt, file.path(OUT_DIR, 'meld_r7b_S4_m5_pairwise.csv'))

# Matrix form for the primary model (first score col alphabetically that includes 'auto',
# else the first one)
prim <- score_cols[grep('auto', score_cols, ignore.case = TRUE)]
if (length(prim) == 0L) prim <- score_cols[[1L]] else prim <- prim[[1L]]
cat(sprintf('\nprimary model for matrix: %s\n', prim))
pm <- pair_dt[model == prim]
mat <- matrix(1, length(LABELS), length(LABELS), dimnames = list(LABELS, LABELS))
for (i in seq_len(nrow(pm))) {
  a <- pm$panel_a[i]; b <- pm$panel_b[i]
  mat[a, b] <- pm$r[i]; mat[b, a] <- pm$r[i]
}
mat_dt <- data.table(panel = LABELS, mat)
fwrite(mat_dt, file.path(OUT_DIR, 'meld_r7b_S4_m5_pairwise_matrix.csv'))

# ---------------------------------------------------------------------------
# Console tables
cat('\n=== S4 M5 intersection summary ===\n')
print(intersect_summary)

cat('\n=== S4 M5 pairwise Pearson r ± 95% block-bootstrap CI (primary model) ===\n')
pm[, txt := sprintf('%.4f [%.4f, %.4f]', r, lo95, hi95)]
print(pm[, .(panel_a, panel_b, txt)])

cat('\n=== S4 M5 correlations to MELD_lambda0 (primary model) ===\n')
ml0 <- pm[panel_a == 'MELD_lambda0' | panel_b == 'MELD_lambda0']
ml0[, other := ifelse(panel_a == 'MELD_lambda0', panel_b, panel_a)]
print(ml0[, .(other, txt)])

cat('\nDONE\n')
