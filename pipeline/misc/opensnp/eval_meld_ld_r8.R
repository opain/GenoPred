#!/usr/bin/env Rscript
# eval_meld_ld_r8.R — Round 8 Section 5 M2 evaluation, chr22.
#
# For one trait, evaluate every candidate LD panel by masked-z re-imputation
# with uniform PSD repair, δ sensitivity, multiple mask draws per block, and
# block bootstrap CIs. Reuses the Round 7b core (m2_core.R).
#
# Candidates:
#   MELD_lambda0_afproj   — analytic target, AF-projection weights (the method)
#   MELD_lambda0_reportN  — same, reported-N weights (near-identical expected)
#   MELD_lambda0_equal    — same, equal weights over the P populations (control)
#   MELD_lambda1_afproj   — pooled target for λ-separation analysis
#   EUR, EAS, AFR, CSA, AMR — single-population empirical (block's R_{POP})
#
# Notes:
#   - MID not evaluated (no block coverage for MID; AF-projection weight for MID
#     is redistributed to EUR in the afproj weight vector).
#   - When a weight vector is missing a population from the candidate set (e.g.
#     reportedN weights for COPD have no CSA), that pop gets w=0 in the target.
#   - Multiple mask draws per block are averaged before block bootstrap CIs.
#
# Usage:  Rscript eval_meld_ld_r8.R <trait>

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
MELD_LD_DIR <- file.path(MISC, 'meld_ld')
source(file.path(MISC, 'm2_core.R'))

args  <- commandArgs(trailingOnly = TRUE)
TRAIT <- args[1]
if (is.na(TRAIT)) stop('usage: eval_meld_ld_r8.R <trait>')

# ---- Config ----------------------------------------------------------------
CHR          <- 22L
POPS         <- c('EUR','EAS','AFR','CSA','AMR')    # block-supported pops
MASK_FRAC    <- 0.10
N_DRAWS_MASK <- 5L
DELTAS       <- c(0.001, 0.01, 0.1)
B_BOOT       <- 2000L
SEED         <- 8000L

# ---- Load harmonised trait data + weights ---------------------------------
d_all <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', TRAIT)))
W_all <- fread(file.path(MISC, 'gbmi_r8_weights.csv'))
W     <- W_all[trait == TRAIT]
if (!nrow(W)) stop(sprintf('no weights for %s', TRAIT))

# Build weight vectors over POPS. MID's AF-projection mass folds into EUR.
mk_w <- function(W_row_source) {
  w <- setNames(rep(0, length(POPS)), POPS)
  for (P in POPS) {
    v <- W[pop == P, get(W_row_source)]
    if (length(v) && is.finite(v)) w[P] <- v
  }
  # Redistribute MID afproj mass into EUR (matches r7 convention).
  if (W_row_source == 'w_afproj') {
    mid <- W[pop == 'MID', w_afproj]
    if (length(mid) && is.finite(mid)) w['EUR'] <- w['EUR'] + mid
  }
  s <- sum(w); if (s > 0) w / s else w
}
w_afproj  <- mk_w('w_afproj')
w_reportN <- mk_w('w_N')
w_equal   <- setNames(rep(1 / length(POPS), length(POPS)), POPS)

cat(sprintf('=== M2 eval: %s (n=%d) ===\n', TRAIT, nrow(d_all)))
cat('afproj weights: ');  print(round(w_afproj, 4))
cat('reportN weights: '); print(round(w_reportN, 4))

# Meta z-scores keyed by rsid.
d_all[, z := meta_beta / meta_se]
z_by_snp <- d_all[, .(SNP = rsid, z)]

# ---- Iterate blocks --------------------------------------------------------
chr_dir <- file.path(MELD_LD_DIR, sprintf('chr%d', CHR))
block_files <- sort(list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE))
stopifnot(length(block_files) > 0)
cat(sprintf('chr%d: %d blocks\n', CHR, length(block_files)))

m2_rows   <- list()
eig_rows  <- list()

for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  if (!all(paste0('R_', POPS) %in% names(bl))) next
  in_ss <- bl$SNP %in% z_by_snp$SNP
  if (sum(in_ss) < 30L) next
  idx <- which(in_ss)
  slice_bl <- list(SNP = bl$SNP[idx])
  for (P in POPS) {
    slice_bl[[paste0('R_', P)]] <- bl[[paste0('R_', P)]][idx, idx, drop = FALSE]
    slice_bl[[paste0('v_', P)]] <- bl[[paste0('v_', P)]][idx]
    slice_bl[[paste0('f_', P)]] <- bl[[paste0('f_', P)]][idx]
  }
  m <- length(idx)

  z_local <- z_by_snp[match(slice_bl$SNP, SNP), z]
  if (any(!is.finite(z_local))) next

  # Build all candidates (pre-repair R, then repair).
  cand_pre <- list()
  cand_pre[['MELD_lambda0_afproj']]  <- reconstruct_R_lambda_P(slice_bl, 0, w_afproj,  POPS)
  cand_pre[['MELD_lambda0_reportN']] <- reconstruct_R_lambda_P(slice_bl, 0, w_reportN, POPS)
  cand_pre[['MELD_lambda0_equal']]   <- reconstruct_R_lambda_P(slice_bl, 0, w_equal,   POPS)
  cand_pre[['MELD_lambda1_afproj']]  <- reconstruct_R_lambda_P(slice_bl, 1, w_afproj,  POPS)
  for (P in POPS) cand_pre[[P]] <- slice_bl[[paste0('R_', P)]]

  cand_R <- lapply(cand_pre, psd_repair)

  # Record eigenvalue diagnostics (min-eig pre/post).
  for (cn in names(cand_R)) {
    eig_rows[[length(eig_rows) + 1L]] <- data.table(
      trait = TRAIT, chr = CHR, block_id = bl$block_id, m = m,
      candidate = cn,
      min_eig_pre  = min_eig(cand_pre[[cn]]),
      min_eig_post = min_eig(cand_R[[cn]])
    )
  }

  # Multiple mask draws.
  set.seed(SEED + fi)
  for (dr in seq_len(N_DRAWS_MASK)) {
    mask_S <- sort(sample.int(m, size = max(2L, round(MASK_FRAC * m))))
    for (delta in DELTAS) {
      for (cn in names(cand_R)) {
        r2 <- run_M2_one(cand_R[[cn]], z_local, mask_S, delta)
        m2_rows[[length(m2_rows) + 1L]] <- data.table(
          trait = TRAIT, chr = CHR, block_id = bl$block_id, m = m,
          candidate = cn, delta = delta, mask_draw = dr, r2 = r2
        )
      }
    }
  }

  if (fi %% 5L == 0L) cat(sprintf('  %d/%d blocks\n', fi, length(block_files)))
}

m2 <- rbindlist(m2_rows)
eig <- rbindlist(eig_rows)

# ---- Average mask draws within block, then block bootstrap ---------------
m2_block <- m2[, .(r2 = mean(r2, na.rm = TRUE), m = first(m)),
               by = .(trait, chr, block_id, candidate, delta)]

set.seed(SEED + 999L)
ci_rows <- list()
for (dl in DELTAS) {
  for (cn in unique(m2_block$candidate)) {
    sub <- m2_block[candidate == cn & delta == dl]
    if (nrow(sub) == 0L) next
    st <- boot_mean(sub$r2, sub$m, B = B_BOOT)
    ci_rows[[length(ci_rows) + 1L]] <- data.table(
      trait = TRAIT, delta = dl, candidate = cn,
      mean = st[['mean']], sd = st[['sd']],
      lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
    )
  }
}
ci <- rbindlist(ci_rows)

# ---- Paired differences: MELD-λ0-afproj vs each single-pop + λ1−λ0 --------
paired_rows <- list()
m2_wide <- dcast(m2_block, trait + chr + block_id + m + delta ~ candidate,
                 value.var = 'r2')
for (dl in DELTAS) {
  sub <- m2_wide[delta == dl]
  base <- 'MELD_lambda0_afproj'
  comparisons <- list()
  for (P in POPS) comparisons[[length(comparisons)+1]] <- c(base, P)
  # equal-weights control
  comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_equal')
  # reportedN vs afproj
  comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_reportN')
  # λ1 vs λ0
  comparisons[[length(comparisons)+1]] <- c('MELD_lambda1_afproj', base)
  for (pn in comparisons) {
    a <- sub[[pn[1]]]; b <- sub[[pn[2]]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5L) next
    st <- boot_paired(a[ok], b[ok], sub$m[ok], B = B_BOOT)
    paired_rows[[length(paired_rows) + 1L]] <- data.table(
      trait = TRAIT, delta = dl,
      comparison = sprintf('%s − %s', pn[1], pn[2]),
      mean_diff = st[['mean_diff']], sd = st[['sd']],
      lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
    )
  }
}
paired <- rbindlist(paired_rows)

# ---- Emit ------------------------------------------------------------------
per_block_file <- file.path(MISC, 'gbmi_r8_m2_per_block.csv')
ci_file        <- file.path(MISC, 'gbmi_r8_m2_ci.csv')
paired_file    <- file.path(MISC, 'gbmi_r8_m2_paired_ci.csv')
eig_file       <- file.path(MISC, 'gbmi_r8_m2_eigen.csv')

# Read existing, drop this trait's rows, append fresh.
append_or_start <- function(dt, fp) {
  if (file.exists(fp)) {
    prior <- fread(fp)
    prior <- prior[trait != TRAIT]
    dt <- rbind(prior, dt, use.names = TRUE, fill = TRUE)
  }
  fwrite(dt, fp)
}
append_or_start(m2_block, per_block_file)
append_or_start(ci,       ci_file)
append_or_start(paired,   paired_file)
append_or_start(eig,      eig_file)

cat(sprintf('\nwrote %s (%d rows)\n', per_block_file, nrow(m2_block)))
cat(sprintf('wrote %s (%d rows)\n',   ci_file,        nrow(ci)))
cat(sprintf('wrote %s (%d rows)\n',   paired_file,    nrow(paired)))
cat(sprintf('wrote %s (%d rows)\n',   eig_file,       nrow(eig)))

# ---- Console summary at δ = 0.01 -----------------------------------------
cat(sprintf('\n=== %s M2 r² (δ = 0.01) ===\n', TRAIT))
sub_ci <- ci[delta == 0.01]
sub_ci[, txt := sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95)]
print(sub_ci[, .(candidate, txt)])

cat(sprintf('\n=== %s paired diffs (δ = 0.01) ===\n', TRAIT))
sub_p <- paired[delta == 0.01]
sub_p[, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(sub_p[, .(comparison, txt)])

cat(sprintf('\n=== %s min-eig summary (post-repair should be ~0) ===\n', TRAIT))
eig_sum <- eig[, .(min_eig_pre_median  = median(min_eig_pre),
                   min_eig_pre_p05     = quantile(min_eig_pre, 0.05),
                   min_eig_post_median = median(min_eig_post)), by = candidate]
print(eig_sum)

cat('\nDONE\n')
