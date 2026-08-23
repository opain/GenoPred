#!/usr/bin/env Rscript
# eval_meld_ld_r9.R — Round 9 Section 4 M2 evaluation on UKB LD panels.
#
# Same protocol as eval_meld_ld_r8.R. Two differences:
#   (a) MELD_LD_DIR points at meld_ld_ukb/chr22 (74 sub-blocks after
#       common-refinement reconciliation);
#   (b) POPS = {EUR, EAS, AFR, CSA, AMR, MID} — MID is now a candidate and
#       is no longer folded into EUR in the afproj weight vector.
#
# Reuses m2_core.R verbatim.
#
# Usage: Rscript eval_meld_ld_r9.R <trait>

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r9')
MELD_LD_DIR <- meld_ld_ukb()
source(file.path(MISC, 'm2_core.R'))

args  <- commandArgs(trailingOnly = TRUE)
TRAIT <- args[1]
if (is.na(TRAIT)) stop('usage: eval_meld_ld_r9.R <trait>')

# ---- Config ----------------------------------------------------------------
CHR          <- 22L
POPS         <- c('EUR','EAS','AFR','CSA','AMR','MID')
MASK_FRAC    <- 0.10
N_DRAWS_MASK <- 5L
DELTAS       <- c(0.001, 0.01, 0.1)
B_BOOT       <- 2000L
SEED         <- 9000L
MIN_BLOCK_M  <- 30L   # matches Round 8 threshold in eval_meld_ld_r8.R

# ---- Load harmonised trait data + Round-8 weights ------------------------
d_all <- readRDS(r8_harmonised(TRAIT))
W_all <- fread(file.path(r_results('r8'), 'gbmi_r8_weights.csv'))
W     <- W_all[trait == TRAIT]
if (!nrow(W)) stop(sprintf('no weights for %s', TRAIT))

# Weight vectors over POPS. MID is NO LONGER folded into EUR.
mk_w <- function(source_col) {
  w <- setNames(rep(0, length(POPS)), POPS)
  for (P in POPS) {
    v <- W[pop == P, get(source_col)]
    if (length(v) && is.finite(v)) w[P] <- v
  }
  s <- sum(w); if (s > 0) w / s else w
}
w_afproj  <- mk_w('w_afproj')
w_reportN <- mk_w('w_N')
w_equal   <- setNames(rep(1 / length(POPS), length(POPS)), POPS)

cat(sprintf('=== R9 M2 eval: %s (n=%d) ===\n', TRAIT, nrow(d_all)))
cat('afproj weights (MID no longer folded):\n'); print(round(w_afproj, 4))
cat('reportedN weights:\n'); print(round(w_reportN, 4))

d_all[, z := meta_beta / meta_se]
z_by_snp <- d_all[, .(SNP = rsid, z)]

# ---- Iterate UKB sub-blocks ------------------------------------------------
chr_dir <- file.path(MELD_LD_DIR, sprintf('chr%d', CHR))
block_files <- sort(list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE))
cat(sprintf('chr%d: %d UKB sub-blocks\n', CHR, length(block_files)))

m2_rows  <- list()
eig_rows <- list()

for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  if (!all(paste0('R_', POPS) %in% names(bl))) next
  in_ss <- bl$SNP %in% z_by_snp$SNP
  if (sum(in_ss) < MIN_BLOCK_M) next
  idx <- which(in_ss)
  slice <- list(SNP = bl$SNP[idx])
  for (P in POPS) {
    slice[[paste0('R_', P)]] <- bl[[paste0('R_', P)]][idx, idx, drop = FALSE]
    slice[[paste0('v_', P)]] <- bl[[paste0('v_', P)]][idx]
    slice[[paste0('f_', P)]] <- bl[[paste0('f_', P)]][idx]
  }
  m <- length(idx)
  z_local <- z_by_snp[match(slice$SNP, SNP), z]
  if (any(!is.finite(z_local))) next

  cand_pre <- list()
  cand_pre[['MELD_lambda0_afproj']]  <- reconstruct_R_lambda_P(slice, 0, w_afproj,  POPS)
  cand_pre[['MELD_lambda0_reportN']] <- reconstruct_R_lambda_P(slice, 0, w_reportN, POPS)
  cand_pre[['MELD_lambda0_equal']]   <- reconstruct_R_lambda_P(slice, 0, w_equal,   POPS)
  cand_pre[['MELD_lambda1_afproj']]  <- reconstruct_R_lambda_P(slice, 1, w_afproj,  POPS)
  for (P in POPS) cand_pre[[P]] <- slice[[paste0('R_', P)]]

  cand_R <- lapply(cand_pre, psd_repair)

  for (cn in names(cand_R)) {
    eig_rows[[length(eig_rows) + 1L]] <- data.table(
      trait = TRAIT, chr = CHR, block_id = bl$block_id, m = m,
      candidate = cn,
      min_eig_pre  = min_eig(cand_pre[[cn]]),
      min_eig_post = min_eig(cand_R[[cn]])
    )
  }

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
  if (fi %% 10L == 0L) cat(sprintf('  %d/%d blocks\n', fi, length(block_files)))
}

m2  <- rbindlist(m2_rows)
eig <- rbindlist(eig_rows)

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

paired_rows <- list()
m2_wide <- dcast(m2_block, trait + chr + block_id + m + delta ~ candidate,
                 value.var = 'r2')
for (dl in DELTAS) {
  sub <- m2_wide[delta == dl]
  base <- 'MELD_lambda0_afproj'
  comparisons <- list()
  for (P in POPS) comparisons[[length(comparisons)+1]] <- c(base, P)
  comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_equal')
  comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_reportN')
  comparisons[[length(comparisons)+1]] <- c('MELD_lambda1_afproj', base)
  # Extra: equal − EUR is the "blending" component per the R9 headline
  comparisons[[length(comparisons)+1]] <- c('MELD_lambda0_equal', 'EUR')
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

per_block_file <- file.path(OUT_DIR, 'gbmi_r9_m2_per_block.csv')
ci_file        <- file.path(OUT_DIR, 'gbmi_r9_m2_ci.csv')
paired_file    <- file.path(OUT_DIR, 'gbmi_r9_m2_paired_ci.csv')
eig_file       <- file.path(OUT_DIR, 'gbmi_r9_m2_eigen.csv')

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

cat(sprintf('\nwrote per_block (%d rows), ci (%d), paired (%d), eigen (%d)\n',
            nrow(m2_block), nrow(ci), nrow(paired), nrow(eig)))

cat(sprintf('\n=== %s R9 M2 r² (δ = 0.01) ===\n', TRAIT))
sub <- ci[delta == 0.01]
sub[, txt := sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95)]
print(sub[, .(candidate, txt)])

cat(sprintf('\n=== %s R9 paired diffs (δ = 0.01) ===\n', TRAIT))
sp <- paired[delta == 0.01]
sp[, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(sp[, .(comparison, txt)])

cat(sprintf('\n=== %s R9 min-eig medians ===\n', TRAIT))
eig_sum <- eig[, .(pre_med  = median(min_eig_pre),
                   pre_min  = min(min_eig_pre),
                   post_med = median(min_eig_post)), by = candidate]
print(eig_sum)
cat('\nDONE\n')
