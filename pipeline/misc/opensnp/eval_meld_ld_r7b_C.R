#!/usr/bin/env Rscript
# Round 7b Section C — stratum-impurity test.
#
# The λ=0 derivation assumes strata are pure populations. Round 3 established
# that the Yengo AFR file carries a ~16% EUR-inferred admixture proportion
# (Round-3 recovery on the `eur00` mixture: p_eur = 0.1562, p_afr = 0.8126,
# p_other = 0.0312). Renormalised to a two-way EUR/AFR mixture:
#   p_eur_in_AFR_stratum = 0.1562 / (0.1562 + 0.8126) = 0.1613  (~= 0.16)
#   p_afr_in_AFR_stratum = 0.8126 / (0.1562 + 0.8126) = 0.8387  (~= 0.84)
# Using **exactly p_eur_in_AFR_stratum = 0.161, p_afr_in_AFR_stratum = 0.839**
# below (reported per prompt instruction).
#
# Method:
#   1. Rebuild R_AFR as a λ=0 mixture of 0.161 EUR + 0.839 AFR (per prompt's
#      literal "λ=0 mixture" wording), i.e. treating the AFR stratum's β itself
#      as an inner IVW meta of two disjoint pure-pop subcohorts. This gives
#      Σ_afr_stratum = 0.161 Cov_EUR + 0.839 Cov_AFR; V analogous; R the
#      normalised form. f_afr_stratum = 0.161 f_EUR + 0.839 f_AFR.
#   2. Reconstruct R* using the ORIGINAL synthetic 50/50 EUR/AFR meta but with
#      the corrected R_AFR (and its v, f) in place of the pure-AFR ingredients.
#   3. Re-run Round 6b D3 (permuted-ΔAF placebo) and D5 (|ΔAF| tertile
#      stratification) on chr22 with uniform PSD repair and block bootstrap.
#   4. Compare the pre-correction residual (Round 6b) with the post-correction
#      residual. If it shrinks: stratum impurity was the mechanism.
#
# Chr22 only, w_eur (outer synthetic) = 0.5 for the D3/D5 headline; also runs
# at 0.25 and 0.75 for completeness.
#
# Outputs (chr22):
#   meld_r7b_C_d3_placebo.csv
#   meld_r7b_C_d5_af_stratified.csv
#   meld_r7b_C_summary_ci.csv       consolidated bootstrap CIs for D3 and D5

suppressPackageStartupMessages({
  library(data.table); library(Matrix)
})

sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
meld_ld_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
out_dir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'

CHR          <- 22L
w_grid       <- c(0.75, 0.50, 0.25)
lambda_grid  <- c(0, 1)
n_perm       <- 8L
n_draws_mask <- 5L                  # Section E
mask_frac    <- 0.10
DELTA        <- 0.01
B_BOOT       <- 2000L
N_TOTAL      <- 500000L

# Round-3 renormalised AFR-stratum composition (see header)
P_EUR_IN_AFR <- 0.161
P_AFR_IN_AFR <- 0.839

# ---------------------------------------------------------------------------
# Correct the AFR ingredients into "AFR stratum" ingredients (λ=0 mixture of
# 0.161 EUR + 0.839 AFR).
correct_afr_stratum <- function(bl) {
  Cov_EUR <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_AFR <- bl$R_AFR * outer(sqrt(bl$v_AFR), sqrt(bl$v_AFR))
  Sigma_AFR_new <- P_EUR_IN_AFR * Cov_EUR + P_AFR_IN_AFR * Cov_AFR
  V_AFR_new     <- P_EUR_IN_AFR * bl$v_EUR + P_AFR_IN_AFR * bl$v_AFR
  inv_sd <- 1 / sqrt(V_AFR_new)
  R_AFR_new <- Sigma_AFR_new * outer(inv_sd, inv_sd)
  R_AFR_new[!is.finite(R_AFR_new)] <- 0
  f_AFR_new <- P_EUR_IN_AFR * bl$f_EUR + P_AFR_IN_AFR * bl$f_AFR
  list(R = R_AFR_new, v = V_AFR_new, f = f_AFR_new)
}

# Reconstruction: outer meta of (EUR, AFR-stratum). afr_ingr is either the
# original (bl$R_AFR, bl$v_AFR, bl$f_AFR) or the corrected list from
# correct_afr_stratum(bl).
reconstruct_R_2pop <- function(bl, lambda, w_eur, afr_R, afr_v, afr_f,
                               N_total = N_TOTAL) {
  w_afr <- 1 - w_eur
  N_eur <- w_eur * N_total; N_afr <- w_afr * N_total
  denom <- N_eur + N_afr
  Cov_EUR <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_AFR <- afr_R    * outer(sqrt(afr_v),   sqrt(afr_v))
  Sigma_star <- (N_eur * Cov_EUR + N_afr * Cov_AFR) / denom
  V_star     <- (N_eur * bl$v_EUR + N_afr * afr_v) / denom
  d <- bl$f_EUR - afr_f
  B_scalar <- 4 * w_eur * w_afr
  Sigma <- Sigma_star + lambda * B_scalar * outer(d, d)
  V     <- V_star     + lambda * B_scalar * d^2
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

# Placebo (D3): permute the outer ΔAF vector and rebuild R at λ=1.
reconstruct_R_2pop_dperm <- function(bl, w_eur, afr_R, afr_v, afr_f, d_perm,
                                     N_total = N_TOTAL) {
  w_afr <- 1 - w_eur
  N_eur <- w_eur * N_total; N_afr <- w_afr * N_total
  denom <- N_eur + N_afr
  Cov_EUR <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_AFR <- afr_R    * outer(sqrt(afr_v),   sqrt(afr_v))
  Sigma_star <- (N_eur * Cov_EUR + N_afr * Cov_AFR) / denom
  V_star     <- (N_eur * bl$v_EUR + N_afr * afr_v) / denom
  B_scalar <- 4 * w_eur * w_afr
  Sigma <- Sigma_star + B_scalar * outer(d_perm, d_perm)
  V     <- V_star     + B_scalar * d_perm^2
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

psd_repair <- function(R) {
  ev <- eigen(R, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  R2 <- ev$vectors %*% (lam * t(ev$vectors))
  d <- sqrt(pmax(diag(R2), 0)); d[d == 0] <- 1
  R2 <- R2 / outer(d, d); diag(R2) <- 1
  R2
}
run_M2_one <- function(R_cand, z, mask_S) {
  M <- setdiff(seq_along(z), mask_S)
  R_MM <- R_cand[M, M, drop = FALSE]
  R_SM <- R_cand[mask_S, M, drop = FALSE]
  z_M <- z[M]; z_S <- z[mask_S]
  ev <- suppressWarnings(eigen(R_MM, symmetric = TRUE))
  U <- ev$vectors; lam <- ev$values
  ridge <- DELTA * abs(sum(lam) / length(lam))
  Utz <- crossprod(U, z_M)
  z_hat <- drop(R_SM %*% (U %*% (Utz / (lam + ridge))))
  r <- suppressWarnings(cor(z_hat, z_S))
  if (is.finite(r)) r^2 else NA_real_
}
mean_over_draws <- function(R_cand, z, mask_list) {
  vals <- sapply(mask_list, function(S) run_M2_one(R_cand, z, S))
  mean(vals, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Load per-w synthetic sumstats
z_by_w <- list()
for (w in w_grid) {
  tag <- sprintf('eur%02d', round(100 * w))
  ss  <- fread(file.path(sumstats_dir, sprintf('yengo_2022_height_mix_%s.txt', tag)))
  ss[, z := beta / standard_error]
  z_by_w[[sprintf('w%02d', round(100 * w))]] <- ss[, .(SNP = variant_id, z)]
}

# ---------------------------------------------------------------------------
# Per-block loop — for each afr_treatment ('original', 'corrected'):
#   D3 rows: chr, block_id, m, w_eur, afr_treatment, perm, r2
#   D5 rows: chr, block_id, m, w_eur, afr_treatment, af_bin, af_bin_lo/hi, r2_l0, r2_l1, diff
chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', CHR))
block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
cat(sprintf('chr%d: %d blocks\n', CHR, length(block_files)))

d3_rows <- list(); d5_rows <- list()
for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  keep <- rep(TRUE, length(bl$SNP))
  for (nm in names(z_by_w)) keep <- keep & bl$SNP %in% z_by_w[[nm]]$SNP
  if (sum(keep) < 30L) next
  idx <- which(keep)
  slice_bl <- list(SNP = bl$SNP[idx],
                   R_EUR = bl$R_EUR[idx, idx, drop = FALSE],
                   R_AFR = bl$R_AFR[idx, idx, drop = FALSE],
                   v_EUR = bl$v_EUR[idx], v_AFR = bl$v_AFR[idx],
                   f_EUR = bl$f_EUR[idx], f_AFR = bl$f_AFR[idx])
  m_use <- length(idx)

  corrected <- correct_afr_stratum(slice_bl)
  treatments <- list(
    original = list(R = slice_bl$R_AFR,   v = slice_bl$v_AFR,   f = slice_bl$f_AFR),
    corrected = list(R = corrected$R,     v = corrected$v,      f = corrected$f)
  )

  z_local <- lapply(names(z_by_w), function(nm) z_by_w[[nm]][match(slice_bl$SNP, SNP), z])
  names(z_local) <- names(z_by_w)
  if (any(sapply(z_local, function(z) any(!is.finite(z))))) next

  set.seed(80000L + fi)
  mask_list <- lapply(seq_len(n_draws_mask), function(d)
    sort(sample.int(m_use, size = max(2L, round(mask_frac * m_use)))))

  for (trt in names(treatments)) {
    tr <- treatments[[trt]]
    for (w in w_grid) {
      wtag <- sprintf('w%02d', round(100 * w))
      z_here <- z_local[[wtag]]
      R_l0_rep <- psd_repair(reconstruct_R_2pop(slice_bl, 0, w, tr$R, tr$v, tr$f))
      R_l1_rep <- psd_repair(reconstruct_R_2pop(slice_bl, 1, w, tr$R, tr$v, tr$f))

      # D3: placebo permutation of outer d
      d_outer <- slice_bl$f_EUR - tr$f
      for (perm_i in seq_len(n_perm)) {
        set.seed(90000L + fi * 100L + perm_i + (trt == 'corrected') * 5000L)
        d_perm <- sample(d_outer)
        R_placebo <- psd_repair(reconstruct_R_2pop_dperm(slice_bl, w, tr$R, tr$v, tr$f, d_perm))
        r2 <- mean_over_draws(R_placebo, z_here, mask_list)
        d3_rows[[length(d3_rows) + 1L]] <- data.table(
          chr = CHR, block_id = bl$block_id, m = m_use,
          afr_treatment = trt, w_eur = w, perm = perm_i, r2 = r2)
      }
      # Also record λ=0 and λ=1 in D3 output for direct comparison
      r2_l0 <- mean_over_draws(R_l0_rep, z_here, mask_list)
      r2_l1 <- mean_over_draws(R_l1_rep, z_here, mask_list)
      d3_rows[[length(d3_rows) + 1L]] <- data.table(
        chr = CHR, block_id = bl$block_id, m = m_use,
        afr_treatment = trt, w_eur = w, perm = 0L, r2 = r2_l0)   # perm=0 sentinel = λ=0
      d3_rows[[length(d3_rows) + 1L]] <- data.table(
        chr = CHR, block_id = bl$block_id, m = m_use,
        afr_treatment = trt, w_eur = w, perm = -1L, r2 = r2_l1)  # perm=-1 sentinel = λ=1

      # D5: |ΔAF| tertile stratification
      abs_d <- abs(d_outer)
      bin_cuts <- quantile(abs_d, probs = seq(0, 1, length.out = 4L))
      bin_id <- as.integer(cut(abs_d, breaks = bin_cuts, include.lowest = TRUE))
      bin_id[is.na(bin_id)] <- 1L
      for (b in seq_len(3L)) {
        idx_b <- which(bin_id == b)
        if (length(idx_b) < 30L) next
        set.seed(100000L + fi * 100L + b + (trt == 'corrected') * 5000L)
        mask_b <- sort(sample.int(length(idx_b),
                                  size = max(2L, round(mask_frac * length(idx_b)))))
        mask_list_b <- lapply(seq_len(n_draws_mask), function(d) {
          set.seed(101000L + fi * 1000L + b * 100L + d)
          sort(sample.int(length(idx_b), size = max(2L, round(mask_frac * length(idx_b)))))
        })
        r2_l0_b <- mean(sapply(mask_list_b, function(S)
          run_M2_one(R_l0_rep[idx_b, idx_b, drop = FALSE], z_here[idx_b], S)), na.rm = TRUE)
        r2_l1_b <- mean(sapply(mask_list_b, function(S)
          run_M2_one(R_l1_rep[idx_b, idx_b, drop = FALSE], z_here[idx_b], S)), na.rm = TRUE)
        d5_rows[[length(d5_rows) + 1L]] <- data.table(
          chr = CHR, block_id = bl$block_id,
          afr_treatment = trt, w_eur = w, af_bin = b, m_bin = length(idx_b),
          af_bin_lo = bin_cuts[b], af_bin_hi = bin_cuts[b + 1L],
          r2_l0 = r2_l0_b, r2_l1 = r2_l1_b, diff = r2_l1_b - r2_l0_b)
      }
    }
  }
  if (fi %% 5L == 0L) cat(sprintf('  %d/%d blocks done\n', fi, length(block_files)))
}

d3 <- rbindlist(d3_rows); d5 <- rbindlist(d5_rows)
fwrite(d3, file.path(out_dir, 'meld_r7b_C_d3_placebo.csv'))
fwrite(d5, file.path(out_dir, 'meld_r7b_C_d5_af_stratified.csv'))

# ---------------------------------------------------------------------------
# Bootstrap CIs
sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
boot_paired <- function(a, b, w, B = B_BOOT) {
  diffs <- a - b; boots <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(length(diffs), length(diffs), replace = TRUE)
    boots[i] <- sw(diffs[idx], w[idx])
  }
  c(mean_diff = sw(diffs, w), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

set.seed(2027L)
# Aggregate D3: per-block λ0, λ1, placebo (mean across perms). Then paired
# comparisons: (λ1 − λ0) and (λ1 − placebo).
d3_bl <- dcast(d3[perm == 0L,   .(chr, block_id, m, afr_treatment, w_eur, r2)],
               chr + block_id + m + afr_treatment + w_eur ~ ., value.var = 'r2')
setnames(d3_bl, '.', 'r2_l0')
d3_l1 <- dcast(d3[perm == -1L,  .(chr, block_id, m, afr_treatment, w_eur, r2)],
               chr + block_id + m + afr_treatment + w_eur ~ ., value.var = 'r2')
setnames(d3_l1, '.', 'r2_l1')
d3_pl <- d3[perm > 0L, .(placebo_r2 = mean(r2, na.rm = TRUE)),
            by = .(chr, block_id, m, afr_treatment, w_eur)]
d3_agg <- Reduce(function(a, b) merge(a, b, by = c('chr','block_id','m','afr_treatment','w_eur')),
                 list(d3_bl, d3_l1, d3_pl))
fwrite(d3_agg, file.path(out_dir, 'meld_r7b_C_d3_aggregated.csv'))

summary_rows <- list()
for (trt in unique(d3_agg$afr_treatment)) {
  for (w in unique(d3_agg$w_eur)) {
    sub <- d3_agg[afr_treatment == trt & w_eur == w]
    if (nrow(sub) == 0L) next
    st1 <- boot_paired(sub$r2_l1, sub$r2_l0, sub$m)
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      afr_treatment = trt, w_eur = w,
      metric = 'MELD-λ1 − MELD-λ0',
      mean_diff = st1[['mean_diff']], sd = st1[['sd']],
      lo95 = st1[['lo95.2.5%']], hi95 = st1[['hi95.97.5%']])
    st2 <- boot_paired(sub$r2_l1, sub$placebo_r2, sub$m)
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      afr_treatment = trt, w_eur = w,
      metric = 'MELD-λ1 − placebo',
      mean_diff = st2[['mean_diff']], sd = st2[['sd']],
      lo95 = st2[['lo95.2.5%']], hi95 = st2[['hi95.97.5%']])
  }
}
# D5: (M2_λ1 − M2_λ0) per bin, aggregated per treatment/w/bin
d5_agg <- d5[, .(diff_mean = sw(diff, m_bin),
                 n_blocks = .N), by = .(afr_treatment, w_eur, af_bin)]
for (trt in unique(d5$afr_treatment)) {
  for (w in unique(d5$w_eur)) {
    for (b in unique(d5$af_bin)) {
      sub <- d5[afr_treatment == trt & w_eur == w & af_bin == b]
      if (nrow(sub) < 5L) next
      # Bootstrap over blocks (diff is per block already)
      diffs <- sub$diff; wts <- sub$m_bin
      boots <- replicate(B_BOOT, {
        idx <- sample.int(length(diffs), length(diffs), replace = TRUE)
        sw(diffs[idx], wts[idx])
      })
      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        afr_treatment = trt, w_eur = w,
        metric = sprintf('D5 tertile %d (M2_λ1 − M2_λ0)', b),
        mean_diff = sw(diffs, wts), sd = sd(boots),
        lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
    }
  }
}
summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, 'meld_r7b_C_summary_ci.csv'))

# ---------------------------------------------------------------------------
cat('\n=== Round-3 AFR-stratum composition used: p_EUR=', P_EUR_IN_AFR,
    ' p_AFR=', P_AFR_IN_AFR, ' ===\n', sep = '')

cat('\n=== D3: MELD-λ1 − MELD-λ0 (paired diff, bootstrap CI) ===\n')
summary_dt[, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(dcast(summary_dt[metric == 'MELD-λ1 − MELD-λ0'],
            w_eur ~ afr_treatment, value.var = 'txt'))

cat('\n=== D3: MELD-λ1 − placebo (paired diff) ===\n')
print(dcast(summary_dt[metric == 'MELD-λ1 − placebo'],
            w_eur ~ afr_treatment, value.var = 'txt'))

cat('\n=== D5: (M2_λ1 − M2_λ0) by |ΔAF| tertile ===\n')
d5_txt <- summary_dt[grepl('^D5 tertile', metric)]
print(dcast(d5_txt, w_eur + afr_treatment ~ metric, value.var = 'txt'))

cat('\nDONE\n')
