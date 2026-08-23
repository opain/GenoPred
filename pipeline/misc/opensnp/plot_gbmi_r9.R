#!/usr/bin/env Rscript
# plot_gbmi_r9.R — Round 9 Section 5/6 figures.
#
# (a) Section 5 two-panel: composition component (MELD-λ0 − equal) and
#     blending component (equal − EUR), each with R8 and R9 bars per trait.
# (b) Section 6.1: per-population M2 improvement (R9 − R8) vs log10 panel-size
#     ratio (UKB N / 1KG+HGDP N).
# (c) Section 6.2: MELD − EUR (R9) vs trait's EUR weight.

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r9')
IMG  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'

# --- Load headline paired data -----------------------------------------------
paired_r9  <- fread(file.path(OUT_DIR, 'gbmi_r9_m2_paired_ci.csv'))
paired_r8b <- fread(file.path(OUT_DIR, 'gbmi_r8b_m2_paired_ci.csv'))

comp_dt <- function(paired_dt, source_lbl) {
  paired_dt[delta == 0.01, .(trait, comparison, mean_diff, lo95, hi95)][
    , source := source_lbl]
}
p_r9  <- comp_dt(paired_r9,  'R9_UKB')
p_r8b <- comp_dt(paired_r8b, 'R8_1KG')

get_pair <- function(paired_dt, cmp_str) {
  paired_dt[grepl(cmp_str, comparison, fixed = TRUE)]
}

# Section 5 two-panel: composition + blending
comp_r9  <- get_pair(p_r9,  'MELD_lambda0_afproj − MELD_lambda0_equal')
comp_r8b <- get_pair(p_r8b, 'MELD_lambda0_afproj − MELD_lambda0_equal')
comp <- rbind(comp_r9, comp_r8b)
comp[, kind := 'composition (MELD-lambda0 minus equal)']

blend_r9  <- get_pair(p_r9,  'MELD_lambda0_equal − EUR')
blend_r8b <- get_pair(p_r8b, 'MELD_lambda0_equal − EUR')
blend <- rbind(blend_r9, blend_r8b)
blend[, kind := 'blending (equal minus EUR)']

both <- rbind(comp, blend)
both$source <- factor(both$source, levels = c('R8_1KG','R9_UKB'))
both$kind   <- factor(both$kind,
                      levels = c('composition (MELD-lambda0 minus equal)',
                                 'blending (equal minus EUR)'))

p5 <- ggplot(both, aes(x = trait, y = mean_diff, fill = source)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                position = position_dodge(width = 0.7), width = 0.25,
                colour = 'grey30') +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  facet_wrap(~ kind, ncol = 2, scales = 'free_y') +
  labs(x = NULL, y = 'M2 r2 difference',
       fill = 'LD source',
       title = 'GBMI chr22: composition vs blending, R8 (1KG+HGDP) vs R9 (UKB)',
       subtitle = 'Same sub-block structure; delta = 0.01; error bars = 95% block-boot CI') +
  scale_fill_manual(values = c('R8_1KG' = '#E69F00', 'R9_UKB' = '#0072B2')) +
  theme_half_open() + background_grid() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

png(file.path(IMG, 'gbmi_r9_composition_vs_blending.png'),
    res = 200, width = 2400, height = 1200, units = 'px')
print(p5)
dev.off()

# --- Section 6.1: per-pop M2 improvement vs log10(UKB N / 1KG N) --------------
ci_r9  <- fread(file.path(OUT_DIR, 'gbmi_r9_m2_ci.csv'))
ci_r8b <- fread(file.path(OUT_DIR, 'gbmi_r8b_m2_ci.csv'))

# Panel size table
ukb_N <- fread(file.path('/users/k1806347/oliverpainfel/Data/ukb/zenodo_14614207',
                          'ukb_panel_sizes.csv'))
setnames(ukb_N, 'N', 'N_UKB')
# 1KG+HGDP N per pop, from meld_context_2 (used across R8 blocks)
onekg_N <- data.table(pop = c('EUR','EAS','AFR','CSA','AMR','MID'),
                      N_1KG = c(665, 737, 688, 675, 412, 136))

Ns <- merge(ukb_N, onekg_N, by = 'pop')
Ns[, ratio := N_UKB / N_1KG]
Ns[, log10_ratio := log10(ratio)]

# Per-pop M2 improvement per trait: R9 minus R8b, same block structure
pops <- c('EUR','EAS','AFR','CSA','AMR')
mm <- ci_r9  [delta == 0.01 & candidate %in% pops, .(trait, candidate, m9 = mean)]
m8 <- ci_r8b [delta == 0.01 & candidate %in% pops, .(trait, candidate, m8 = mean)]
per_pop <- merge(m8, mm, by = c('trait','candidate'))
per_pop[, delta := m9 - m8]
setnames(per_pop, 'candidate', 'pop')
per_pop <- merge(per_pop, Ns, by = 'pop')

# Mean improvement per pop, across traits
agg <- per_pop[, .(mean_delta = mean(delta),
                   se = sd(delta) / sqrt(.N),
                   log10_ratio = mean(log10_ratio)), by = pop]

p61 <- ggplot(agg, aes(x = log10_ratio, y = mean_delta, label = pop)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = mean_delta - 1.96*se, ymax = mean_delta + 1.96*se),
                width = 0.05, colour = 'grey30') +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 4.5, force = 3) +
  labs(x = 'log10 panel-size ratio (UKB / 1KG+HGDP)',
       y = 'Per-population M2 improvement (R9 − R8)  averaged over traits',
       title = 'Panel-size boost per population',
       subtitle = 'Larger UKB panels improve small-panel populations more than EUR') +
  theme_half_open() + background_grid()
png(file.path(IMG, 'gbmi_r9_panel_size_boost.png'),
    res = 200, width = 2000, height = 1400, units = 'px')
print(p61)
dev.off()

# --- Section 6.2: MELD − EUR (R9) vs trait's EUR weight ----------------------
paired_diff_EUR_R9 <- paired_r9[delta == 0.01 &
                                comparison == 'MELD_lambda0_afproj − EUR',
                                .(trait, adv_R9 = mean_diff, lo95, hi95)]
# EUR weights (afproj)
W_r8 <- fread(file.path(r_results('r8'), 'gbmi_r8_weights.csv'))
eur_w <- W_r8[pop == 'EUR', .(trait, eur_w = w_afproj)]
d62 <- merge(paired_diff_EUR_R9, eur_w, by = 'trait')

# Also show R8b for context
paired_diff_EUR_R8b <- paired_r8b[delta == 0.01 &
                                  comparison == 'MELD_lambda0_afproj − EUR',
                                  .(trait, adv_R8 = mean_diff)]
d62 <- merge(d62, paired_diff_EUR_R8b, by = 'trait')

p62 <- ggplot(d62, aes(x = eur_w)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.006, colour = 'grey30') +
  geom_point(aes(y = adv_R9), size = 4, colour = '#0072B2') +
  geom_point(aes(y = adv_R8), size = 3, colour = '#E69F00', shape = 17) +
  ggrepel::geom_text_repel(aes(y = adv_R9, label = trait),
                           size = 3.5, force = 2) +
  labs(x = "Trait's EUR weight  (afproj)",
       y = 'M2(MELD-lambda0) minus M2(EUR)',
       title = 'MELD advantage over EUR vs EUR dominance',
       subtitle = 'blue circles = R9 (UKB), orange triangles = R8 (1KG+HGDP), same sub-block structure') +
  theme_half_open() + background_grid()
png(file.path(IMG, 'gbmi_r9_adv_vs_eur_weight.png'),
    res = 200, width = 2000, height = 1400, units = 'px')
print(p62)
dev.off()

cat('=== R9 headline figures written ===\n')
cat('  gbmi_r9_composition_vs_blending.png\n')
cat('  gbmi_r9_panel_size_boost.png\n')
cat('  gbmi_r9_adv_vs_eur_weight.png\n')
