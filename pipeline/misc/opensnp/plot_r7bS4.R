#!/usr/bin/env Rscript
# Round 7b Section 4 figures for the write-up.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(cowplot)
})

OUT <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
FIG <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Target r² per panel (with 95% CI)
marg <- fread(file.path(OUT, 'meld_r7b_S4_target_r2_marginal.csv'))
lvl  <- c('MELD_lambda0','EUR','AFR','AMR','CSA','EAS')
marg[, panel_f := factor(panel, levels = lvl)]
pal <- c(MELD_lambda0 = '#B2182B', EUR = '#2166AC', AFR = '#4393C3',
         AMR = '#92C5DE', CSA = '#F4A582', EAS = '#D6604D')

p1 <- ggplot(marg, aes(x = panel_f, y = r, fill = panel_f, colour = panel_f)) +
  geom_col(alpha = 0.85, width = 0.65) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.2, colour = 'black') +
  geom_text(aes(label = sprintf('%.3f', r), y = hi95 + 0.02), size = 3, colour = 'black') +
  scale_fill_manual(values = pal, guide = 'none') +
  scale_colour_manual(values = pal, guide = 'none') +
  labs(x = NULL, y = 'Pearson r  ( height ~ PGS,  95% bootstrap CI )',
       title = 'Round 7b S4 — OpenSNP EUR height prediction by LD panel',
       subtitle = 'LDpred2-auto on Yengo 2022 all-ancestry, n = 777') +
  theme_half_open() + background_grid()

png(file.path(FIG, 'meld_r7bS4_target_r2.png'),
    res = 200, width = 1800, height = 1200, units = 'px')
print(p1); dev.off()

# ---------------------------------------------------------------------------
# Paired Δr vs MELD-λ0
paired <- fread(file.path(OUT, 'meld_r7b_S4_target_r2_paired_diff.csv'))
paired[, other := sub('^MELD_lambda0 . ', '', comparison)]
paired[, other := sub('.*[−-] ', '', comparison)]
paired[, other_f := factor(other, levels = c('EUR','AFR','AMR','CSA','EAS'))]

p2 <- ggplot(paired, aes(x = other_f, y = diff_r)) +
  geom_hline(yintercept = 0, colour = 'grey40', linetype = 2) +
  geom_col(fill = '#B2182B', alpha = 0.85, width = 0.55) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.2, colour = 'black') +
  geom_text(aes(label = sprintf('%+.3f', diff_r), y = hi95 + 0.005),
            size = 3, colour = 'black') +
  labs(x = NULL, y = expression(Delta*r == r[MELD-lambda0] - r[panel]),
       title = 'Round 7b S4 — paired Δr vs MELD-λ0 (same-individual bootstrap)') +
  theme_half_open() + background_grid()

png(file.path(FIG, 'meld_r7bS4_target_r2_paired.png'),
    res = 200, width = 1800, height = 1000, units = 'px')
print(p2); dev.off()

# ---------------------------------------------------------------------------
# Cross-panel PGS correlation matrix as a tile plot
mat_dt <- fread(file.path(OUT, 'meld_r7b_S4_target_pgs_matrix.csv'))
m <- as.matrix(mat_dt[, -'panel'])
rownames(m) <- mat_dt$panel
long <- as.data.table(as.table(m))
setnames(long, c('a','b','r'))
long[, a := factor(a, levels = lvl)]
long[, b := factor(b, levels = lvl)]

p3 <- ggplot(long, aes(x = a, y = b, fill = r)) +
  geom_tile(colour = 'white') +
  geom_text(aes(label = sprintf('%.3f', r)), size = 3.5) +
  scale_fill_gradient(low = 'white', high = '#B2182B', limits = c(0.6, 1),
                      name = 'r') +
  labs(x = NULL, y = NULL,
       title = 'Round 7b S4 — pairwise Pearson r of per-individual PGS across panels',
       subtitle = 'OpenSNP EUR (n = 777), LDpred2-auto') +
  theme_minimal_grid()

png(file.path(FIG, 'meld_r7bS4_pgs_correlation_matrix.png'),
    res = 200, width = 1600, height = 1400, units = 'px')
print(p3); dev.off()

cat('wrote:\n')
cat('  ', file.path(FIG, 'meld_r7bS4_target_r2.png'), '\n')
cat('  ', file.path(FIG, 'meld_r7bS4_target_r2_paired.png'), '\n')
cat('  ', file.path(FIG, 'meld_r7bS4_pgs_correlation_matrix.png'), '\n')
