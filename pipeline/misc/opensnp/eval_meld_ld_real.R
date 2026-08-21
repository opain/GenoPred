#!/usr/bin/env Rscript
# Round 7 Section 2 — real-sumstats M2 on Yengo 2022 height all-ancestry meta.
#
# Extends the two-population MELD-λ0 target to P populations. Candidates on M2:
#   MELD_lambda0     analytic P-population meta target R*
#   EUR, EAS, AFR, CSA, AMR   single-population empirical R (chr22 has all five)
#   MELD_lambda1     P-population mega target (rank-(P-1) B), scored once to close it out
#
# Weights are computed two ways per the prompt:
#   (a) from per-ancestry median N in the Yengo per-pop files
#   (b) from bigsnpr::snp_ancestry_summary AF-projection on Yengo-all EAF
#
# Round 6b fixes are applied: uniform PSD repair (eigen-clamp + diag rescale)
# to every candidate before scoring, block-bootstrap 95% CIs on the M2 mean.
# No M1 (no ground truth on real data).
#
# Chr22 only. Outputs:
#   meld_r7_s2_weights.csv                weight vectors (a) and (b) and discrepancy
#   meld_r7_s2_m2_per_block.csv           per-block r² per candidate per weight-source
#   meld_r7_s2_m2_ci.csv                  bootstrap-CI M2 summary
#   meld_r7_s2_m2_paired_ci.csv           paired candidate differences with CIs

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
  library(bigreadr)
})
.libPaths(c('/home/claude/Rlibs', .libPaths()))

# ---------------------------------------------------------------------------
sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
meld_ld_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
emp_var_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
BIGSNPR_DIR  <- '/users/k1806347/oliverpainfel/Data/bigsnpr'
out_dir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'

CHR          <- 22L                # chr22 only (see prompt: full-genome deferred to Section 5)
pops         <- c('EUR', 'EAS', 'AFR', 'CSA', 'AMR')
mask_frac    <- 0.10
n_draws_mask <- 1L
DELTA        <- 0.01
B_BOOT       <- 2000L

# ---------------------------------------------------------------------------
# 1. Per-ancestry N weights from Yengo per-pop files
cat('=== weights (a): per-ancestry median N from Yengo per-pop files ===\n')
per_pop_files <- c(EUR = 'yengo_2022_height_eur.txt',
                   EAS = 'yengo_2022_height_eas.txt',
                   AFR = 'yengo_2022_height_afr.txt',
                   CSA = 'yengo_2022_height_sas.txt',    # SAS ↔ CSA per meld_context
                   AMR = 'yengo_2022_height_amr.txt')
N_by_pop <- sapply(per_pop_files, function(fn) {
  d <- fread(file.path(sumstats_dir, fn), select = 'n')
  as.numeric(median(d$n))
})
w_paperN <- N_by_pop / sum(N_by_pop)
cat('per-pop median N:\n'); print(N_by_pop)
cat('weights from paper N:\n'); print(round(w_paperN, 4))

# ---------------------------------------------------------------------------
# 2. AF-projection weights (bigsnpr::snp_ancestry_summary) on Yengo-all EAF
cat('\n=== weights (b): AF-projection recovery from Yengo-all EAF ===\n')
all_freq   <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'ref_freqs.csv.gz'))
projection <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'projection.csv.gz'))
# Privé's ancestry vignette correction (constant)
correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

coarse_group <- function(fine) {
  g <- fine
  g[g %in% c('Scandinavia','United Kingdom','Ireland')]     <- 'Europe (North West)'
  g[g %in% c('Europe (South East)','Europe (North East)')]  <- 'Europe (East)'
  g
}
fine_pops <- colnames(all_freq)[-(1:5)]
grp_fct   <- factor(coarse_group(fine_pops), levels = unique(coarse_group(fine_pops)))
super_map <- c(
  'Africa (West)'='AFR','Africa (South)'='AFR','Africa (East)'='AFR','Africa (North)'='AFR',
  'Middle East'='MID','Ashkenazi'='EUR','Italy'='EUR','Finland'='EUR',
  'Europe (East)'='EUR','Europe (North West)'='EUR','Europe (South West)'='EUR',
  'South America'='AMR','Sri Lanka'='CSA','Pakistan'='CSA','Bangladesh'='CSA',
  'Asia (East)'='EAS','Japan'='EAS','Philippines'='EAS'
)

ss_all <- fread(file.path(sumstats_dir, 'yengo_2022_height_all.txt'))
gwas_freq <- ss_all[, .(chr = as.integer(chromosome),
                        pos = base_pair_location,
                        rsid = variant_id, a0 = other_allele, a1 = effect_allele,
                        freq = effect_allele_frequency)]
gwas_freq[, beta := 1]                 # snp_match needs a beta col
matched <- snp_match(as.data.frame(gwas_freq), all_freq[, 1:5], match.min.prop = 0.05)
matched$freq <- ifelse(matched$beta < 0, 1 - matched$freq, matched$freq)

res <- snp_ancestry_summary(
  freq          = matched$freq,
  info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
  projection    = projection[matched$`_NUM_ID_`, -(1:5)],
  correction    = correction
)
by_coarse <- tapply(res, grp_fct, sum)
by_super  <- tapply(by_coarse, super_map[names(by_coarse)], sum, default = 0)
w_afproj  <- sapply(c('EUR','EAS','AFR','CSA','AMR','MID'), function(s)
                    if (s %in% names(by_super)) as.numeric(by_super[[s]]) else 0)
# Renormalise dropping MID (we don't have a MID single-pop candidate)
if ('MID' %in% names(w_afproj) && w_afproj[['MID']] > 0) {
  cat(sprintf('AF-proj: MID recovered at %.4f — reassigning to nearest (EUR) so weights match candidate set\n', w_afproj[['MID']]))
  w_afproj[['EUR']] <- w_afproj[['EUR']] + w_afproj[['MID']]
  w_afproj <- w_afproj[names(w_afproj) != 'MID']
}
w_afproj <- w_afproj[pops]
w_afproj <- w_afproj / sum(w_afproj)
cat('weights from AF-projection:\n'); print(round(w_afproj, 4))

weights_dt <- data.table(
  pop           = pops,
  N_median      = N_by_pop[pops],
  w_paperN      = w_paperN[pops],
  w_afproj      = w_afproj[pops],
  discrepancy   = w_afproj[pops] - w_paperN[pops]
)
fwrite(weights_dt, file.path(out_dir, 'meld_r7_s2_weights.csv'))
cat(sprintf('wrote %s\n', file.path(out_dir, 'meld_r7_s2_weights.csv')))

# ---------------------------------------------------------------------------
# 3. Reconstruction helper (P populations)
reconstruct_R_lambda_P <- function(bl, lambda, w_pop, pops_here) {
  # Weighted meta covariance
  Sigma_star <- matrix(0, length(bl$SNP), length(bl$SNP))
  V_star     <- rep(0, length(bl$SNP))
  for (p in pops_here) {
    R_p  <- bl[[paste0('R_', p)]]
    v_p  <- bl[[paste0('v_', p)]]
    Cov_p <- R_p * outer(sqrt(v_p), sqrt(v_p))
    Sigma_star <- Sigma_star + w_pop[[p]] * Cov_p
    V_star     <- V_star     + w_pop[[p]] * v_p
  }
  if (abs(lambda) > 1e-12) {
    # B for P populations: Σ_p w_p (f_p - f_meta)(f_p - f_meta)^T at each SNP pair, scaled.
    # Wahlund form: B(i,j) = 4 Σ_p w_p (f_p,i - f_bar_i)(f_p,j - f_bar_j)
    # where f_bar = Σ_p w_p f_p. Also on-diagonal: B(i,i) = 4 Σ_p w_p (f_p,i - f_bar_i)^2.
    F <- do.call(cbind, lapply(pops_here, function(p) bl[[paste0('f_', p)]]))
    fbar <- as.vector(F %*% w_pop[pops_here])
    dF <- F - fbar
    W  <- diag(w_pop[pops_here])
    # B = 4 * dF %*% W %*% t(dF). Rank <= P.
    B <- 4 * (dF %*% W %*% t(dF))
    Sigma <- Sigma_star + lambda * B
    V     <- V_star     + lambda * diag(B)
  } else {
    Sigma <- Sigma_star
    V     <- V_star
  }
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
  z_M  <- z[M]; z_S <- z[mask_S]
  ev <- suppressWarnings(eigen(R_MM, symmetric = TRUE))
  U <- ev$vectors; lam <- ev$values
  ridge <- DELTA * abs(sum(lam) / length(lam))
  Utz <- crossprod(U, z_M)
  z_hat <- drop(R_SM %*% (U %*% (Utz / (lam + ridge))))
  r <- suppressWarnings(cor(z_hat, z_S))
  if (is.finite(r)) r^2 else NA_real_
}

# ---------------------------------------------------------------------------
# 4. Load block RDS for chr22 and iterate
chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', CHR))
block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
cat(sprintf('\nchr%d: %d blocks\n', CHR, length(block_files)))

# Yengo-all z-scores (rebuilt from beta/SE)
ss_all[, z := beta / standard_error]
z_by_snp <- ss_all[, .(SNP = variant_id, z)]

# Two weight sets to test
weight_sets <- list(paperN = w_paperN[pops], afproj = w_afproj[pops])

m2_rows <- list()
for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  # Sanity: must have all 5 pops
  if (!all(paste0('R_', pops) %in% names(bl))) {
    cat(sprintf('  skip block %d — missing pops (chr%d likely not extended)\n', bl$block_id, bl$chr))
    next
  }

  # Restrict to SNPs also in Yengo-all
  in_ss <- bl$SNP %in% z_by_snp$SNP
  if (sum(in_ss) < 30L) next
  idx <- which(in_ss)
  slice_bl <- list(SNP = bl$SNP[idx])
  for (p in pops) {
    slice_bl[[paste0('R_', p)]] <- bl[[paste0('R_', p)]][idx, idx, drop = FALSE]
    slice_bl[[paste0('v_', p)]] <- bl[[paste0('v_', p)]][idx]
    slice_bl[[paste0('f_', p)]] <- bl[[paste0('f_', p)]][idx]
  }
  m_use <- length(idx)

  z_local <- z_by_snp[match(slice_bl$SNP, SNP), z]
  if (any(!is.finite(z_local))) next

  # Mask (shared across candidates and weight-sets for identifiability)
  set.seed(7000L + fi)
  mask_S <- sort(sample.int(m_use, size = max(2L, round(mask_frac * m_use))))

  for (ws in names(weight_sets)) {
    w_here <- weight_sets[[ws]]
    R_star <- psd_repair(reconstruct_R_lambda_P(slice_bl, 0, w_here, pops))
    R_pool <- psd_repair(reconstruct_R_lambda_P(slice_bl, 1, w_here, pops))

    cand_R <- list(MELD_lambda0 = R_star, MELD_lambda1 = R_pool)
    for (p in pops) cand_R[[p]] <- psd_repair(slice_bl[[paste0('R_', p)]])

    for (cn in names(cand_R)) {
      r2 <- run_M2_one(cand_R[[cn]], z_local, mask_S)
      m2_rows[[length(m2_rows) + 1L]] <- data.table(
        chr = CHR, block_id = bl$block_id, m = m_use,
        weight_set = ws, candidate = cn, r2 = r2
      )
    }
  }

  if (fi %% 5L == 0L) cat(sprintf('  %d/%d blocks done\n', fi, length(block_files)))
}

m2 <- rbindlist(m2_rows)
fwrite(m2, file.path(out_dir, 'meld_r7_s2_m2_per_block.csv'))
cat(sprintf('\nwrote %s (%d rows)\n',
            file.path(out_dir, 'meld_r7_s2_m2_per_block.csv'), nrow(m2)))

# ---------------------------------------------------------------------------
# 5. Bootstrap CI on M2 and paired candidate differences
sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
boot_mean <- function(vals, weights, B = B_BOOT) {
  n <- length(vals)
  boots <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    boots[b] <- sw(vals[idx], weights[idx])
  }
  c(mean = sw(vals, weights), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}
boot_paired <- function(a, b, w, B = B_BOOT) {
  diffs <- a - b
  boots <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(length(diffs), length(diffs), replace = TRUE)
    boots[i] <- sw(diffs[idx], w[idx])
  }
  c(mean_diff = sw(diffs, w), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

set.seed(2027L)
ci_rows <- list()
for (ws in unique(m2$weight_set)) {
  for (cn in unique(m2$candidate)) {
    sub <- m2[weight_set == ws & candidate == cn]
    if (nrow(sub) == 0L) next
    st <- boot_mean(sub$r2, sub$m)
    ci_rows[[length(ci_rows) + 1L]] <- data.table(
      weight_set = ws, candidate = cn,
      mean = st[['mean']], sd = st[['sd']],
      lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
    )
  }
}
ci_dt <- rbindlist(ci_rows)
fwrite(ci_dt, file.path(out_dir, 'meld_r7_s2_m2_ci.csv'))

# Paired: MELD-λ0 − each single-pop, and MELD-λ1 − MELD-λ0
paired_rows <- list()
m2_wide <- dcast(m2, chr + block_id + m + weight_set ~ candidate, value.var = 'r2')
for (ws in unique(m2_wide$weight_set)) {
  sub <- m2_wide[weight_set == ws]
  for (cn in setdiff(names(cand_R <- list()), NULL)) NULL
  comparisons <- list(
    c('MELD_lambda0', 'EUR'),
    c('MELD_lambda0', 'EAS'),
    c('MELD_lambda0', 'AFR'),
    c('MELD_lambda0', 'CSA'),
    c('MELD_lambda0', 'AMR'),
    c('MELD_lambda1', 'MELD_lambda0')
  )
  for (pn in comparisons) {
    a <- sub[[pn[1]]]; b <- sub[[pn[2]]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5L) next
    st <- boot_paired(a[ok], b[ok], sub$m[ok])
    paired_rows[[length(paired_rows) + 1L]] <- data.table(
      weight_set = ws, comparison = sprintf('%s − %s', pn[1], pn[2]),
      mean_diff = st[['mean_diff']], sd = st[['sd']],
      lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
    )
  }
}
paired_dt <- rbindlist(paired_rows)
fwrite(paired_dt, file.path(out_dir, 'meld_r7_s2_m2_paired_ci.csv'))

# ---------------------------------------------------------------------------
# Console summaries
cat('\n=== S2 weights ===\n'); print(weights_dt)

cat('\n=== S2 M2 r² by candidate × weight-set (mean [95% CI], post-repair, δ=0.01) ===\n')
ci_dt[, txt := sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95)]
print(dcast(ci_dt, candidate ~ weight_set, value.var = 'txt'))

cat('\n=== S2 paired differences (95% CI) ===\n')
paired_dt[, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(dcast(paired_dt, comparison ~ weight_set, value.var = 'txt'))

cat('\nDONE\n')
