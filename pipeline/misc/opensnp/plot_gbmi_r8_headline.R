#!/usr/bin/env Rscript
# plot_gbmi_r8_headline.R — Round 8 Section 6.
#
# Two headline scatter plots on chr22:
#   §6.1  M2(MELD-λ0-afproj) − M2(best single-pop) vs effective ancestry groups
#         with bootstrap CIs. Yengo height added as an anchor at (1.6, +0.027)
#         if the Round 7 CSV is available.
#   §6.2  M2(MELD-λ1-afproj) − M2(MELD-λ0-afproj) vs mean between-pop variance.
#
# Inputs:  gbmi_r8_m2_ci.csv, gbmi_r8_m2_paired_ci.csv, gbmi_r8_composition.csv,
#          meld_r7_s2_m2_ci.csv (optional, Yengo anchor)
# Outputs: docs/Images/OpenSNP/gbmi_r8_headline_composition.png
#          docs/Images/OpenSNP/gbmi_r8_headline_lambda.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

MISC   <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
IMG    <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
DELTA  <- 0.01

ci     <- fread(file.path(MISC, 'gbmi_r8_m2_ci.csv'))
paired <- fread(file.path(MISC, 'gbmi_r8_m2_paired_ci.csv'))
comp   <- fread(file.path(MISC, 'gbmi_r8_composition.csv'))

# ---- §6.1  M2(MELD-λ0) − M2(best single-pop) vs effective ancestry groups ----
best_sp <- ci[delta == DELTA & candidate %in% c('EUR','EAS','AFR','CSA','AMR'),
              .(best_sp = max(mean), best_sp_lo = max(lo95), best_sp_hi = max(hi95),
                best_sp_name = candidate[which.max(mean)]),
              by = trait]
meld0   <- ci[delta == DELTA & candidate == 'MELD_lambda0_afproj',
              .(trait, meld0 = mean, meld0_lo = lo95, meld0_hi = hi95)]

# For the paired difference with CI, use the paired-CI file rather than
# recomputing from unpaired CIs (the paired CIs are tighter).
paired_meld0 <- paired[delta == DELTA]
paired_meld0[, comp_from := trimws(sapply(strsplit(comparison, '−'), `[`, 1))]
paired_meld0[, comp_to   := trimws(sapply(strsplit(comparison, '−'), `[`, 2))]
# Get M2(λ0-afproj) − each single pop, per trait
diffs_sp <- paired_meld0[comp_from == 'MELD_lambda0_afproj' &
                          comp_to %in% c('EUR','EAS','AFR','CSA','AMR'),
                          .(trait, pop = comp_to, diff = mean_diff,
                            lo95, hi95)]
# For each trait: pick the *smallest* diff (i.e., the single-pop closest to MELD-λ0
# from above; MELD-λ0 − best single-pop is the smallest positive diff).
diff_best <- diffs_sp[, .SD[which.min(diff)], by = trait]
setnames(diff_best, 'pop', 'best_sp')

d61 <- merge(diff_best, comp[, .(trait, eff_groups, mean_B_diag)], by = 'trait')
setnames(d61, c('diff','lo95','hi95'), c('adv','adv_lo','adv_hi'))
fwrite(d61, file.path(MISC, 'gbmi_r8_headline_61.csv'))

p61 <- ggplot(d61, aes(x = eff_groups, y = adv)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = adv_lo, ymax = adv_hi), width = 0.03, colour = 'grey30') +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = paste0(trait, ' (', best_sp, ')')),
                           size = 3, max.overlaps = 20, force = 2) +
  labs(x = 'Effective number of ancestry groups   1 / sum(w_p^2)',
       y = 'M2(MELD lambda0)  minus  M2(best single-pop panel)',
       title = 'GBMI chr22: composition benefit vs balance',
       subtitle = sprintf('delta = %g; error bars = 95%% block-bootstrap CI', DELTA)) +
  theme_half_open() + background_grid()

# Optional: add Yengo height anchor (approx 1.6 eff groups, +0.027 advantage)
yengo_r7 <- file.path(MISC, 'meld_r7_s2_m2_ci.csv')
if (file.exists(yengo_r7)) {
  y7 <- fread(yengo_r7)
  # Round 7 file has weight_set 'afproj' and candidate names EUR/EAS/AFR/CSA/AMR/MELD_lambda0
  y_meld <- y7[weight_set == 'afproj' & candidate == 'MELD_lambda0', mean][1]
  y_sp   <- y7[weight_set == 'afproj' & candidate %in% c('EUR','EAS','AFR','CSA','AMR'), max(mean)]
  y_adv  <- y_meld - y_sp
  # Approx eff_groups from R7 known composition (0.78/0.13/0.05/0.02/0.02) = 1.6
  y_eff  <- 1 / sum(c(0.78, 0.127, 0.053, 0.019, 0.021)^2)
  p61 <- p61 +
    annotate('point', x = y_eff, y = y_adv, shape = 17, colour = 'firebrick', size = 3.5) +
    annotate('text', x = y_eff, y = y_adv, label = 'Yengo height',
             hjust = -0.15, vjust = 0, colour = 'firebrick', size = 3)
}

png(file.path(IMG, 'gbmi_r8_headline_composition.png'),
    res = 200, width = 2000, height = 1400, units = 'px')
print(p61)
dev.off()

# ---- §6.2  M2(λ1) − M2(λ0)  vs  mean between-pop variance -------------------
paired_lam <- paired[delta == DELTA & comparison == 'MELD_lambda1_afproj − MELD_lambda0_afproj',
                     .(trait, mean_diff, lo95, hi95)]
d62 <- merge(paired_lam, comp[, .(trait, eff_groups, mean_B_diag)], by = 'trait')
fwrite(d62, file.path(MISC, 'gbmi_r8_headline_62.csv'))

p62 <- ggplot(d62, aes(x = mean_B_diag, y = mean_diff)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.001, colour = 'grey30') +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = trait), size = 3, max.overlaps = 20, force = 2) +
  labs(x = 'Mean between-population variance   mean_i 4 sum_p w_p (f_p,i - fbar_i)^2',
       y = 'M2(MELD lambda1)  minus  M2(MELD lambda0)',
       title = 'GBMI chr22: lambda separation vs between-population variance',
       subtitle = sprintf('delta = %g; error bars = 95%% block-bootstrap CI', DELTA)) +
  theme_half_open() + background_grid()

png(file.path(IMG, 'gbmi_r8_headline_lambda.png'),
    res = 200, width = 2000, height = 1400, units = 'px')
print(p62)
dev.off()

cat('=== headline data — §6.1 ===\n')
print(d61)
cat('\n=== headline data — §6.2 ===\n')
print(d62)
cat('DONE\n')
