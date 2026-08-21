#!/usr/bin/env Rscript
# Round 7b Section A — weight sweep (A2) + wrong-weight controls (A1) on real
# Yengo-2022 height all-ancestry sumstats, chr22, with the Round 6b conditioning
# fix and Section E's stronger noise handling (5 mask draws averaged per block
# before block-bootstrap).
#
# Purpose: separate the "correct composition" effect from the generic "any
# mixture regularises noisy matrices" effect that the +0.027 r² MELD-λ0 vs EUR
# margin in Round 7 cannot itself distinguish.
#
# Outputs (chr22):
#   meld_r7b_A_weights_used.csv       weight vectors evaluated
#   meld_r7b_A_m2_per_block.csv       per-block per-vector r² (averaged over draws)
#   meld_r7b_A_m2_ci.csv              bootstrap CIs per weight vector
#   meld_r7b_A_paired_ci.csv          paired diffs vs EUR and vs correct MELD-λ0
#   docs/Images/OpenSNP/meld_r7b_sweep.png

suppressPackageStartupMessages({
  library(data.table); library(bigsnpr); library(bigstatsr); library(Matrix)
  library(bigreadr);   library(ggplot2); library(cowplot)
})
.libPaths(c('/home/claude/Rlibs', .libPaths()))

sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
meld_ld_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
BIGSNPR_DIR  <- '/users/k1806347/oliverpainfel/Data/bigsnpr'
out_dir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
figs_dir     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

CHR          <- 22L
pops         <- c('EUR', 'EAS', 'AFR', 'CSA', 'AMR')
mask_frac    <- 0.10
n_draws_mask <- 5L                  # Section E: average draws within-block
DELTA        <- 0.01
B_BOOT       <- 2000L
sweep_grid   <- seq(0, 1, by = 0.1)

# ---------------------------------------------------------------------------
# 1. Recompute the two ground-truth weight vectors, reporting MID pre-fold
per_pop_files <- c(EUR = 'yengo_2022_height_eur.txt',
                   EAS = 'yengo_2022_height_eas.txt',
                   AFR = 'yengo_2022_height_afr.txt',
                   CSA = 'yengo_2022_height_sas.txt',
                   AMR = 'yengo_2022_height_amr.txt')
N_by_pop <- sapply(per_pop_files, function(fn) {
  as.numeric(median(fread(file.path(sumstats_dir, fn), select = 'n')$n))
})
w_paperN <- N_by_pop / sum(N_by_pop)

all_freq   <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'ref_freqs.csv.gz'))
projection <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'projection.csv.gz'))
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
                        freq = effect_allele_frequency, beta = 1)]
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
w_afproj_pre_fold <- sapply(c('EUR','EAS','AFR','CSA','AMR','MID'), function(s)
                            if (s %in% names(by_super)) as.numeric(by_super[[s]]) else 0)
cat(sprintf('AF-projection weights BEFORE MID fold:\n'))
print(round(w_afproj_pre_fold, 4))
w_mid <- w_afproj_pre_fold[['MID']]
cat(sprintf('MID weight (folded into EUR): %.4f\n', w_mid))

w_afproj <- w_afproj_pre_fold
w_afproj[['EUR']] <- w_afproj[['EUR']] + w_afproj[['MID']]
w_afproj <- w_afproj[names(w_afproj) != 'MID']
w_afproj <- w_afproj / sum(w_afproj)

# ---------------------------------------------------------------------------
# 2. Assemble weight vectors to evaluate

# Base relative shares among non-EUR pops from afproj
non_eur <- setdiff(pops, 'EUR')
non_eur_share <- w_afproj[non_eur] / sum(w_afproj[non_eur])
mk_sweep_w <- function(w_eur) {
  w <- c(EUR = w_eur, non_eur_share * (1 - w_eur))
  w[pops]
}
sweep_weights <- lapply(sweep_grid, mk_sweep_w)
names(sweep_weights) <- sprintf('sweep_wEUR_%.2f', sweep_grid)

wrong_weights <- list(
  equal = setNames(rep(0.2, 5), pops),
  # Scramble 1: swap EUR ↔ AFR (dominant → minor)
  scramble_swap_EUR_AFR = {
    w <- w_afproj; e <- w[['EUR']]; w[['EUR']] <- w[['AFR']]; w[['AFR']] <- e; w
  },
  # Scramble 2: cyclic permutation EUR → EAS → AFR → CSA → AMR → EUR
  scramble_cyclic = {
    w <- w_afproj
    v <- c(w[['AMR']], w[['EUR']], w[['EAS']], w[['AFR']], w[['CSA']])
    setNames(v, pops)
  }
)

true_weights <- list(
  paperN = w_paperN[pops],
  afproj = w_afproj[pops]
)

all_weights <- c(true_weights, wrong_weights, sweep_weights)
weights_used_dt <- rbindlist(lapply(names(all_weights), function(nm)
  data.table(label = nm, pop = pops, w = as.numeric(all_weights[[nm]][pops]))))
fwrite(weights_used_dt, file.path(out_dir, 'meld_r7b_A_weights_used.csv'))
cat('\nAll weight vectors evaluated (long-form):\n')
print(dcast(weights_used_dt, label ~ pop, value.var = 'w')[, .(label, EUR, EAS, AFR, CSA, AMR)])

# ---------------------------------------------------------------------------
# 3. Reconstruction / M2 helpers (P populations)
reconstruct_R_star_P <- function(bl, w_pop) {
  Sigma_star <- matrix(0, length(bl$SNP), length(bl$SNP))
  V_star     <- rep(0, length(bl$SNP))
  for (p in pops) {
    R_p  <- bl[[paste0('R_', p)]]
    v_p  <- bl[[paste0('v_', p)]]
    Cov_p <- R_p * outer(sqrt(v_p), sqrt(v_p))
    Sigma_star <- Sigma_star + w_pop[[p]] * Cov_p
    V_star     <- V_star     + w_pop[[p]] * v_p
  }
  inv_sd <- 1 / sqrt(V_star)
  R <- Sigma_star * outer(inv_sd, inv_sd)
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
# 4. Per-block loop
chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', CHR))
block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
cat(sprintf('\nchr%d: %d blocks\n', CHR, length(block_files)))
ss_all[, z := beta / standard_error]
z_by_snp <- ss_all[, .(SNP = variant_id, z)]

rows <- list()
for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  if (!all(paste0('R_', pops) %in% names(bl))) next
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

  # Section E: 5 mask draws averaged within-block (same masks across all candidates)
  set.seed(70000L + fi)
  mask_list <- lapply(seq_len(n_draws_mask), function(d)
    sort(sample.int(m_use, size = max(2L, round(mask_frac * m_use)))))

  # Also score EUR-alone (weight-independent) as a fixed baseline
  R_EUR_rep <- psd_repair(slice_bl$R_EUR)
  r2_EUR <- mean_over_draws(R_EUR_rep, z_local, mask_list)
  rows[[length(rows) + 1L]] <- data.table(
    chr = CHR, block_id = bl$block_id, m = m_use,
    label = 'EUR_baseline', kind = 'single_pop',
    r2 = r2_EUR)

  # Every weight vector => MELD-λ0 reconstruction and M2
  for (nm in names(all_weights)) {
    w_here <- all_weights[[nm]]
    R_star_rep <- psd_repair(reconstruct_R_star_P(slice_bl, w_here))
    r2 <- mean_over_draws(R_star_rep, z_local, mask_list)
    kind <- if (nm %in% names(true_weights))  'true_weights' else
            if (nm %in% names(wrong_weights)) 'wrong_weight' else 'sweep'
    rows[[length(rows) + 1L]] <- data.table(
      chr = CHR, block_id = bl$block_id, m = m_use,
      label = nm, kind = kind, r2 = r2)
  }

  if (fi %% 5L == 0L) cat(sprintf('  %d/%d blocks done\n', fi, length(block_files)))
}

m2 <- rbindlist(rows)
fwrite(m2, file.path(out_dir, 'meld_r7b_A_m2_per_block.csv'))
cat(sprintf('\nwrote %s (%d rows)\n',
            file.path(out_dir, 'meld_r7b_A_m2_per_block.csv'), nrow(m2)))

# ---------------------------------------------------------------------------
# 5. Bootstrap CIs
sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}
boot_mean <- function(vals, weights, B = B_BOOT) {
  n <- length(vals); boots <- numeric(B)
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
for (nm in unique(m2$label)) {
  sub <- m2[label == nm]
  if (nrow(sub) == 0L) next
  st <- boot_mean(sub$r2, sub$m)
  ci_rows[[length(ci_rows) + 1L]] <- data.table(
    label = nm, kind = sub$kind[1L],
    mean = st[['mean']], sd = st[['sd']],
    lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
  )
}
ci_dt <- rbindlist(ci_rows)
fwrite(ci_dt, file.path(out_dir, 'meld_r7b_A_m2_ci.csv'))

# Paired diffs vs EUR baseline and vs afproj MELD-λ0 (the "correct" weight)
paired_rows <- list()
m2_wide <- dcast(m2, chr + block_id + m ~ label, value.var = 'r2')
correct_lbl <- 'afproj'
if (!(correct_lbl %in% names(m2_wide))) correct_lbl <- 'sweep_wEUR_0.80'  # nearest fallback

for (nm in unique(m2$label)) {
  if (nm == 'EUR_baseline') next
  a <- m2_wide[[nm]]; b <- m2_wide[['EUR_baseline']]
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5L) next
  st <- boot_paired(a[ok], b[ok], m2_wide$m[ok])
  paired_rows[[length(paired_rows) + 1L]] <- data.table(
    label = nm, vs = 'EUR_baseline',
    mean_diff = st[['mean_diff']], sd = st[['sd']],
    lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']])

  if (nm != correct_lbl) {
    a2 <- m2_wide[[correct_lbl]]
    ok2 <- is.finite(a) & is.finite(a2)
    if (sum(ok2) >= 5L) {
      st2 <- boot_paired(a[ok2], a2[ok2], m2_wide$m[ok2])
      paired_rows[[length(paired_rows) + 1L]] <- data.table(
        label = nm, vs = correct_lbl,
        mean_diff = st2[['mean_diff']], sd = st2[['sd']],
        lo95 = st2[['lo95.2.5%']], hi95 = st2[['hi95.97.5%']])
    }
  }
}
paired_dt <- rbindlist(paired_rows)
fwrite(paired_dt, file.path(out_dir, 'meld_r7b_A_paired_ci.csv'))

# ---------------------------------------------------------------------------
# 6. Sweep plot with CIs
sweep_ci <- ci_dt[kind == 'sweep']
sweep_ci[, w_EUR := as.numeric(sub('sweep_wEUR_', '', label))]
setorder(sweep_ci, w_EUR)
eur_ci   <- ci_dt[label == 'EUR_baseline']
afproj_ci <- ci_dt[label == 'afproj']
paperN_ci <- ci_dt[label == 'paperN']

p <- ggplot(sweep_ci, aes(x = w_EUR, y = mean)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = '#DDDDDD') +
  geom_line(size = 1) + geom_point(size = 2.5) +
  geom_hline(yintercept = eur_ci$mean, colour = '#2166AC', linetype = 'dashed') +
  geom_hline(yintercept = eur_ci$lo95, colour = '#2166AC', linetype = 'dotted', alpha = .5) +
  geom_hline(yintercept = eur_ci$hi95, colour = '#2166AC', linetype = 'dotted', alpha = .5) +
  geom_vline(xintercept = 0.7800, colour = '#B2182B', linetype = 'solid', alpha = .7) +
  geom_vline(xintercept = 0.7723, colour = '#B2182B', linetype = 'dotted', alpha = .7) +
  annotate('text', x = 0.7800, y = min(sweep_ci$lo95), hjust = -0.05, vjust = 0,
           label = 'w_EUR = 0.78 (afproj)', colour = '#B2182B', size = 3) +
  annotate('text', x = 0.7723, y = max(sweep_ci$hi95), hjust = 1.05, vjust = 1,
           label = 'w_EUR = 0.77 (paper N)', colour = '#B2182B', size = 3) +
  annotate('text', x = 0.05, y = eur_ci$mean, hjust = 0, vjust = -0.5,
           label = sprintf('EUR baseline = %.3f', eur_ci$mean), colour = '#2166AC', size = 3) +
  labs(x = 'w_EUR (remaining rescaled from AF-projection shares)',
       y = 'M2 masked-z r² (post-repair, δ=0.01, size-weighted mean, 95% bootstrap CI)',
       title = 'Round 7b A2: MELD-λ0 M2 weight sweep on real Yengo-2022 height',
       subtitle = sprintf('chr22, %d blocks, %d mask draws averaged per block',
                          length(unique(m2$block_id)), n_draws_mask)) +
  theme_half_open() + background_grid()

png(file.path(figs_dir, 'meld_r7b_sweep.png'),
    width = 2000, height = 1400, res = 200)
print(p); dev.off()

# ---------------------------------------------------------------------------
# Console summaries
cat('\n=== Weight vectors reported ===\n')
print(round(t(sapply(all_weights, function(w) w[pops])), 4))

cat('\n=== MID weight recovered by AF-projection (pre-fold): ',
    round(w_mid, 4), ' — folded into EUR for the candidate set\n')

cat('\n=== M2 CIs by weight vector (post-repair, δ=0.01, averaged over ',
    n_draws_mask, ' draws) ===\n', sep = '')
ci_dt[, txt := sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95)]
print(ci_dt[order(-mean), .(label, kind, txt)])

cat('\n=== Paired vs EUR baseline ===\n')
paired_dt[vs == 'EUR_baseline', txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(paired_dt[vs == 'EUR_baseline', .(label, txt)])

cat('\n=== Paired vs correct-weight MELD-λ0 (', correct_lbl, ') ===\n', sep = '')
paired_dt[vs == correct_lbl, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(paired_dt[vs == correct_lbl, .(label, txt)])

cat('\nDONE\n')
