#!/usr/bin/env Rscript
# Round 6 step 3 checkpoint figures: M1 (Frobenius error per block) and
# M2 (masked-z re-imputation r²) per candidate, per mixture point.
#
# Produces:
#   docs/Images/OpenSNP/meld_r6_m1_hist.png
#   docs/Images/OpenSNP/meld_r6_m1_by_w.png
#   docs/Images/OpenSNP/meld_r6_m2_by_w.png
#
# Prints per-candidate summary tables to stdout for the check-back to the user.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

opensnp_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
figs_dir    <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

m1 <- fread(file.path(opensnp_dir, 'meld_r6_m1.csv'))
m2 <- fread(file.path(opensnp_dir, 'meld_r6_m2.csv'))

cand_order <- c('MELD_lambda0', 'MELD_lambda1', 'MELD_pooled',
                'EUR', 'AFR', 'EAS', 'CSA', 'AMR')
m1[, candidate := factor(candidate, levels = cand_order)]
m2[, candidate := factor(candidate, levels = cand_order)]

# Palette: MELD variants in reds, single-pop in cool colours
pal <- c(MELD_lambda0 = '#B2182B', MELD_lambda1 = '#EF8A62', MELD_pooled = '#FDDBC7',
         EUR = '#2166AC', AFR = '#4393C3', EAS = '#92C5DE',
         CSA = '#8073AC', AMR = '#B2ABD2')

# ---------------------------------------------------------------------------
# M1 histogram at w = 0.5
p1 <- ggplot(m1[w_eur == 0.5], aes(x = err, fill = candidate)) +
  geom_histogram(bins = 60, alpha = 0.85, colour = NA) +
  facet_wrap(~ candidate, scales = 'free_y') +
  scale_fill_manual(values = pal) +
  labs(x = expression('||R'[cand]*' - R*||'[F]*' / ||R*||'[F]),
       y = 'blocks',
       title = 'M1: block-wise Frobenius error to R* at w_eur = 0.50') +
  theme_half_open() + background_grid() + theme(legend.position = 'none')
png(file.path(figs_dir, 'meld_r6_m1_hist.png'),
    width = 2000, height = 1500, res = 200)
print(p1); dev.off()

# ---------------------------------------------------------------------------
# M1 mean vs w by candidate
m1_sum <- m1[, .(err_mean = mean(err, na.rm = TRUE),
                 err_se   = sd(err,  na.rm = TRUE) / sqrt(sum(!is.na(err)))),
             by = .(candidate, w_eur)]
p2 <- ggplot(m1_sum, aes(x = factor(w_eur), y = err_mean,
                          colour = candidate, group = candidate)) +
  geom_line(size = 1) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = err_mean - err_se, ymax = err_mean + err_se), width = 0.15) +
  scale_colour_manual(values = pal) +
  labs(x = 'w_eur', y = 'mean Frobenius err to R*',
       title = 'M1: mean block-wise error by mixture point') +
  theme_half_open() + background_grid()
png(file.path(figs_dir, 'meld_r6_m1_by_w.png'),
    width = 2000, height = 1200, res = 200)
print(p2); dev.off()

# ---------------------------------------------------------------------------
# M2 mean r² vs w by candidate (per canonical δ = 0.01)
m2_sum <- m2[delta == 0.01,
             .(r2_mean_wm = sum(r2 * m, na.rm = TRUE) / sum(m * (!is.na(r2))),
               r2_mean    = mean(r2, na.rm = TRUE)),
             by = .(candidate, w_eur)]
p3 <- ggplot(m2_sum, aes(x = factor(w_eur), y = r2_mean_wm,
                         colour = candidate, group = candidate)) +
  geom_line(size = 1) + geom_point(size = 3) +
  scale_colour_manual(values = pal) +
  labs(x = 'w_eur', y = 'M2: mean masked-z r² (size-weighted, δ = 0.01)',
       title = 'M2: mean masked-z re-imputation r² by mixture point') +
  theme_half_open() + background_grid()
png(file.path(figs_dir, 'meld_r6_m2_by_w.png'),
    width = 2000, height = 1200, res = 200)
print(p3); dev.off()

# ---------------------------------------------------------------------------
# Print summary tables
cat('\n=== M1 summary (mean err by candidate × w) ===\n')
print(dcast(m1_sum, candidate ~ w_eur, value.var = 'err_mean'))
cat('\n=== M2 summary at δ=0.01 (size-weighted mean r² by candidate × w) ===\n')
print(dcast(m2_sum, candidate ~ w_eur, value.var = 'r2_mean_wm'))
cat('\nDONE\n')
