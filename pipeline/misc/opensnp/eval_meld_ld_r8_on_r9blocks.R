#!/usr/bin/env Rscript
# eval_meld_ld_r8_on_r9blocks.R — R8 M2 evaluated on the R9 sub-block
# structure, for the R8-vs-R9 paired comparison at fixed block structure.
#
# Same M2 protocol as eval_meld_ld_r8.R, but MELD_LD_DIR points at
# meld_ld_r8_on_r9blocks/chr22 (67 sub-blocks with R8's 1KG+HGDP LD data
# restricted to the R9 sub-block SNP set). POPS excludes MID since R8 has
# no MID data.

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r9')
MELD_LD_DIR <- meld_ld_r8_on_r9blocks()
source(file.path(MISC, 'm2_core.R'))

args  <- commandArgs(trailingOnly = TRUE)
TRAIT <- args[1]
if (is.na(TRAIT)) stop('usage: eval_meld_ld_r8_on_r9blocks.R <trait>')

CHR          <- 22L
POPS         <- c('EUR','EAS','AFR','CSA','AMR')
MASK_FRAC    <- 0.10
N_DRAWS_MASK <- 5L
DELTAS       <- c(0.001, 0.01, 0.1)
B_BOOT       <- 2000L
SEED         <- 8500L
MIN_BLOCK_M  <- 30L

d_all <- readRDS(r8_harmonised(TRAIT))
W_all <- fread(file.path(r_results('r8'), 'gbmi_r8_weights.csv'))
W     <- W_all[trait == TRAIT]
if (!nrow(W)) stop(sprintf('no weights for %s', TRAIT))

mk_w <- function(source_col) {
  w <- setNames(rep(0, length(POPS)), POPS)
  for (P in POPS) {
    v <- W[pop == P, get(source_col)]
    if (length(v) && is.finite(v)) w[P] <- v
  }
  # For R8-style: MID mass folds into EUR (matches R8 convention).
  if (source_col == 'w_afproj') {
    mid <- W[pop == 'MID', w_afproj]
    if (length(mid) && is.finite(mid)) w['EUR'] <- w['EUR'] + mid
  }
  s <- sum(w); if (s > 0) w / s else w
}
w_afproj  <- mk_w('w_afproj')
w_reportN <- mk_w('w_N')
w_equal   <- setNames(rep(1 / length(POPS), length(POPS)), POPS)

cat(sprintf('=== R8-on-R9blocks: %s (n=%d) ===\n', TRAIT, nrow(d_all)))
d_all[, z := meta_beta / meta_se]
z_by_snp <- d_all[, .(SNP = rsid, z)]

chr_dir <- file.path(MELD_LD_DIR, sprintf('chr%d', CHR))
block_files <- sort(list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE))
cat(sprintf('chr%d: %d R8-on-R9 blocks\n', CHR, length(block_files)))

m2_rows  <- list()
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
}
m2 <- rbindlist(m2_rows)
m2_block <- m2[, .(r2 = mean(r2, na.rm = TRUE), m = first(m)),
               by = .(trait, chr, block_id, candidate, delta)]

set.seed(SEED + 999L)
ci_rows <- list()
for (dl in DELTAS) for (cn in unique(m2_block$candidate)) {
  sub <- m2_block[candidate == cn & delta == dl]
  if (nrow(sub) == 0L) next
  st <- boot_mean(sub$r2, sub$m, B = B_BOOT)
  ci_rows[[length(ci_rows) + 1L]] <- data.table(
    trait = TRAIT, delta = dl, candidate = cn,
    mean = st[['mean']], sd = st[['sd']],
    lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
  )
}
ci <- rbindlist(ci_rows)

paired_rows <- list()
m2_wide <- dcast(m2_block, trait + chr + block_id + m + delta ~ candidate, value.var = 'r2')
for (dl in DELTAS) {
  sub <- m2_wide[delta == dl]
  base <- 'MELD_lambda0_afproj'
  comps <- list()
  for (P in POPS) comps[[length(comps)+1]] <- c(base, P)
  comps[[length(comps)+1]] <- c(base, 'MELD_lambda0_equal')
  comps[[length(comps)+1]] <- c(base, 'MELD_lambda0_reportN')
  comps[[length(comps)+1]] <- c('MELD_lambda1_afproj', base)
  comps[[length(comps)+1]] <- c('MELD_lambda0_equal', 'EUR')
  for (pn in comps) {
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

append_or_start <- function(dt, fp) {
  if (file.exists(fp)) {
    prior <- fread(fp)
    prior <- prior[trait != TRAIT]
    dt <- rbind(prior, dt, use.names = TRUE, fill = TRUE)
  }
  fwrite(dt, fp)
}
append_or_start(m2_block, file.path(OUT_DIR, 'gbmi_r8b_m2_per_block.csv'))
append_or_start(ci,       file.path(OUT_DIR, 'gbmi_r8b_m2_ci.csv'))
append_or_start(paired,   file.path(OUT_DIR, 'gbmi_r8b_m2_paired_ci.csv'))

cat(sprintf('DONE %s: per_block=%d, ci=%d, paired=%d\n',
            TRAIT, nrow(m2_block), nrow(ci), nrow(paired)))
