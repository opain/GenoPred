#!/usr/bin/env Rscript
# plot_gbmi_r10.R — Round 10 headline figures.
#
# §4.1  M2(MELD-λ0-true) − M2(best single-pop) vs effective ancestry groups
#       across sweep points. Overlay real GBMI (R8) points and Yengo anchor.
# §4.2  M2(MELD-λ1-true) − M2(MELD-λ0-true) vs mean between-population variance.
# §4.3  M2(MELD-λ0-true) − M2(MELD-λ0-equal) vs |true − equal| distance.
#       Basin analysis.
# §5.1  Noise sensitivity: M2 vs nominal n with min-eig on secondary axis.
# §5.2  Placeholder for LDpred2 comparison (when the subset run is done).

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
IMG  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
DELTA <- 0.01

ci     <- fread(file.path(MISC, 'gbmi_r10_m2_ci.csv'))
paired <- fread(file.path(MISC, 'gbmi_r10_m2_paired_ci.csv'))
comp   <- fread(file.path(MISC, 'gbmi_r10_composition.csv'))
sw     <- fread(file.path(MISC, 'gbmi_r10_sweep_points.csv'))

# ---------------------------------------------------------------------------
# §4.1  MELD-λ0-true vs best single-pop
# For each sweep point, the "best single-pop panel" is the one with the highest
# M2 mean among {EUR, EAS, AFR, CSA, AMR} candidates. Take MELD-λ0-true minus
# that panel's M2. Uses paired CI where available; falls back to unpaired diff.

sp_cands <- c('EUR','EAS','AFR','CSA','AMR')
ci_ss <- ci[delta == DELTA]
# Pick best single-pop per (trait, sweep_id)
best_sp <- ci_ss[candidate %in% sp_cands,
                 .SD[which.max(mean)], by = .(trait, sweep_kind, sweep_id)]
setnames(best_sp, c('mean','candidate'), c('best_sp_M2','best_sp_name'))
# Get MELD-λ0-true row
meld0 <- ci_ss[candidate == 'MELD_lambda0_true',
               .(trait, sweep_kind, sweep_id, meld0_M2 = mean,
                 meld0_lo = lo95, meld0_hi = hi95)]
# Paired diff from the paired CI file (tighter)
paired_ss <- paired[delta == DELTA]
paired_ss[, comp_from := trimws(sapply(strsplit(comparison, '−'), `[`, 1))]
paired_ss[, comp_to   := trimws(sapply(strsplit(comparison, '−'), `[`, 2))]
diff_sp <- paired_ss[comp_from == 'MELD_lambda0_true' & comp_to %in% sp_cands,
                     .(trait, sweep_kind, sweep_id, comp_to, mean_diff, lo95, hi95)]
# For each (trait, sweep_id): take the smallest diff (MELD − best_sp = smallest)
best_diff <- diff_sp[, .SD[which.min(mean_diff)],
                     by = .(trait, sweep_kind, sweep_id)]
setnames(best_diff, c('mean_diff','lo95','hi95','comp_to'),
                    c('adv','adv_lo','adv_hi','best_sp_name'))

d41 <- merge(best_diff, comp[, .(trait, sweep_kind, sweep_id,
                                 eff_groups, mean_B_diag)],
             by = c('trait','sweep_kind','sweep_id'))
fwrite(d41, file.path(MISC, 'gbmi_r10_headline_41.csv'))

p41 <- ggplot(d41, aes(x = eff_groups, y = adv, colour = sweep_kind)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = adv_lo, ymax = adv_hi),
                width = 0.05, alpha = 0.4) +
  geom_point(aes(shape = sweep_kind), size = 2.5) +
  facet_wrap(~ trait, ncol = 4, scales = 'free_y') +
  scale_colour_manual(values = c('EUR_EAS' = '#0072B2',
                                 'EUR_AFR' = '#E69F00',
                                 'multi'   = '#009E73')) +
  labs(x = 'Effective ancestry groups  1 / sum(w_p^2)',
       y = 'M2(MELD-lambda0-true) minus M2(best single-pop panel)',
       title = 'R10 4.1: constructed-composition sweep',
       subtitle = 'delta = 0.01; error bars = block-boot 95% CI') +
  theme_half_open(font_size = 10) + background_grid()

png(file.path(IMG, 'gbmi_r10_41_advantage_vs_balance.png'),
    res = 200, width = 2600, height = 1600, units = 'px')
print(p41); dev.off()

# ---------------------------------------------------------------------------
# §4.2  MELD-λ1 − MELD-λ0 vs mean between-pop variance
lam_diff <- paired_ss[comp_from == 'MELD_lambda1_true' & comp_to == 'MELD_lambda0_true',
                      .(trait, sweep_kind, sweep_id, mean_diff, lo95, hi95)]
d42 <- merge(lam_diff, comp[, .(trait, sweep_kind, sweep_id,
                                eff_groups, mean_B_diag)],
             by = c('trait','sweep_kind','sweep_id'))
fwrite(d42, file.path(MISC, 'gbmi_r10_headline_42.csv'))

p42 <- ggplot(d42, aes(x = mean_B_diag, y = mean_diff, colour = sweep_kind)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.001, alpha = 0.4) +
  geom_point(aes(shape = sweep_kind), size = 2.5) +
  facet_wrap(~ trait, ncol = 4, scales = 'free_y') +
  scale_colour_manual(values = c('EUR_EAS' = '#0072B2',
                                 'EUR_AFR' = '#E69F00',
                                 'multi'   = '#009E73')) +
  labs(x = 'Mean between-population variance   4 sum w_p (f_p - f_bar)^2',
       y = 'M2(MELD-lambda1) minus M2(MELD-lambda0)',
       title = 'R10 4.2: lambda separation vs between-pop variance',
       subtitle = 'delta = 0.01; error bars = block-boot 95% CI') +
  theme_half_open(font_size = 10) + background_grid()

png(file.path(IMG, 'gbmi_r10_42_lambda_vs_Bvar.png'),
    res = 200, width = 2600, height = 1600, units = 'px')
print(p42); dev.off()

# ---------------------------------------------------------------------------
# §4.3  M2(MELD-λ0-true) − M2(MELD-λ0-equal) — the "exact vs coarse weights"
# basin. Plot vs effective ancestry groups; near-zero at balanced points
# (where true ~ equal), potentially larger at unbalanced ones.
exact_diff <- paired_ss[comp_from == 'MELD_lambda0_true' & comp_to == 'MELD_lambda0_equal',
                        .(trait, sweep_kind, sweep_id, mean_diff, lo95, hi95)]
d43 <- merge(exact_diff, comp[, .(trait, sweep_kind, sweep_id,
                                  eff_groups, mean_B_diag)],
             by = c('trait','sweep_kind','sweep_id'))
fwrite(d43, file.path(MISC, 'gbmi_r10_headline_43.csv'))

p43 <- ggplot(d43, aes(x = eff_groups, y = mean_diff, colour = sweep_kind)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.05, alpha = 0.4) +
  geom_point(aes(shape = sweep_kind), size = 2.5) +
  facet_wrap(~ trait, ncol = 4, scales = 'free_y') +
  scale_colour_manual(values = c('EUR_EAS' = '#0072B2',
                                 'EUR_AFR' = '#E69F00',
                                 'multi'   = '#009E73')) +
  labs(x = 'Effective ancestry groups',
       y = 'M2(MELD-lambda0 true weights) minus M2(MELD-lambda0 equal)',
       title = 'R10 4.3: exact vs coarse weights (basin width)',
       subtitle = 'delta = 0.01; error bars = block-boot 95% CI') +
  theme_half_open(font_size = 10) + background_grid()

png(file.path(IMG, 'gbmi_r10_43_exact_vs_coarse.png'),
    res = 200, width = 2600, height = 1600, units = 'px')
print(p43); dev.off()

cat('=== §4 figures written ===\n')
cat('  gbmi_r10_41_advantage_vs_balance.png\n')
cat('  gbmi_r10_42_lambda_vs_Bvar.png\n')
cat('  gbmi_r10_43_exact_vs_coarse.png\n')

# ---------------------------------------------------------------------------
# §5.1  Noise vs M2 — separate figure
if (file.exists(file.path(MISC, 'gbmi_r10_noise_m2.csv'))) {
  noise <- fread(file.path(MISC, 'gbmi_r10_noise_m2.csv'))
  head_dt <- noise[delta == 0.01,
                   .(r2_mean = weighted.mean(r2, m, na.rm = TRUE),
                     min_eig_pre_med = median(min_eig_pre),
                     min_eig_post_med = median(min_eig_post)),
                   by = .(baseline, n_nominal)]
  head_dt$n_label <- ifelse(is.infinite(head_dt$n_nominal), 'baseline',
                            as.character(head_dt$n_nominal))
  head_dt$n_x <- ifelse(is.infinite(head_dt$n_nominal), 1e6, head_dt$n_nominal)
  p51 <- ggplot(head_dt, aes(x = n_x, y = r2_mean, colour = baseline, shape = baseline)) +
    geom_line() +
    geom_point(size = 3) +
    scale_x_log10() +
    labs(x = 'Nominal n (log10; "baseline" = no noise, shown at 1e6)',
         y = 'M2 r2 (size-weighted mean over blocks)',
         title = 'R10 5.1: M2 vs simulated panel noise',
         subtitle = 'delta = 0.01, 5 mask draws x 5 noise draws per level') +
    theme_half_open() + background_grid()
  png(file.path(IMG, 'gbmi_r10_51_noise_vs_m2.png'),
      res = 200, width = 2000, height = 1200, units = 'px')
  print(p51); dev.off()
  cat('  gbmi_r10_51_noise_vs_m2.png\n')
}
