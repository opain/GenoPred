#!/usr/bin/env Rscript
# gbmi_r10_sanity.R — Round 10 Section 1 gate.
#
# Reconstruct each trait's constructed z at N_tilde = real reported per-arm
# median N and compare against the real GBMI multi-ancestry meta z-scores.
# Per prompt §1: "It will not be exact — GBMI meta-analyses per biobank
# then across ancestries, and applies its own QC — but it should be high.
# If it is not, the construction is wrong; stop and report."
#
# Reports:
#   - Pearson r(z_synth, z_meta) per trait
#   - Distribution of differences |z_synth - z_meta|
#   - Coverage: fraction of harmonised SNPs where both z's are finite
#
# Emits: gbmi_r10_sanity.csv, docs/Images/OpenSNP/gbmi_r10_sanity.png

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
IMG  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/docs/Images/OpenSNP'
source(file.path(MISC, 'gbmi_r10_synth.R'))

TRAITS <- c('Asthma','COPD','Gout','HF','IPF','Stroke','VTE')
ARM_POPS <- list(
  Asthma = c('EUR','EAS','AFR','CSA','AMR'),
  COPD   = c('EUR','EAS','AFR','AMR'),
  Gout   = c('EUR','EAS','AFR','AMR'),
  HF     = c('EUR','EAS','AFR','AMR'),
  IPF    = c('EUR','EAS','AFR','AMR'),
  Stroke = c('EUR','EAS','AFR','AMR'),
  VTE    = c('EUR','EAS','AFR','AMR')
)
THRESHOLD <- 0.95  # per plan

rows <- list()
plots <- list()

for (TR in TRAITS) {
  d <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', TR)))
  pops <- ARM_POPS[[TR]]
  # Real reported per-arm median N (used as N_tilde for the sanity check)
  N_rep <- r10_reported_N(d, pops)

  synth <- r10_synth(d, N_rep)
  # Real meta z from the harmonised RDS
  z_meta <- d$meta_beta / d$meta_se
  ok <- synth$ok & is.finite(z_meta) & is.finite(synth$z_synth)
  n_ok <- sum(ok)
  r <- cor(synth$z_synth[ok], z_meta[ok])
  diff <- synth$z_synth[ok] - z_meta[ok]

  cat(sprintf('=== %s ===\n', TR))
  cat(sprintf('  n_ok = %d / %d\n', n_ok, nrow(d)))
  cat(sprintf('  r(z_synth, z_meta) = %.4f\n', r))
  cat(sprintf('  median |diff| = %.4f, p95 = %.4f, max = %.4f\n',
              median(abs(diff)), quantile(abs(diff), 0.95), max(abs(diff))))
  cat('  reported N per arm: '); print(round(N_rep))

  rows[[length(rows) + 1L]] <- data.table(
    trait = TR, n_snps = nrow(d), n_ok = n_ok,
    r_zsynth_zmeta = r,
    median_abs_diff = median(abs(diff)),
    p95_abs_diff = quantile(abs(diff), 0.95),
    max_abs_diff = max(abs(diff)),
    N_reported = paste(sprintf('%s=%d', pops, round(N_rep)), collapse = ';')
  )

  # Scatter for the plot panel
  set.seed(10001L)
  n_pt <- min(n_ok, 2000L)
  ix <- sample(which(ok), n_pt)
  plots[[TR]] <- data.table(trait = TR,
                            z_meta = z_meta[ix],
                            z_synth = synth$z_synth[ix])
}

res <- rbindlist(rows)
fwrite(res, file.path(MISC, 'gbmi_r10_sanity.csv'))
cat(sprintf('\nwrote %s\n', file.path(MISC, 'gbmi_r10_sanity.csv')))

# Scatter facets
pl <- rbindlist(plots)
p <- ggplot(pl, aes(x = z_meta, y = z_synth)) +
  geom_abline(slope = 1, intercept = 0, colour = 'firebrick', linetype = 'dashed') +
  geom_point(alpha = 0.25, size = 0.6) +
  facet_wrap(~ trait, ncol = 4, scales = 'free') +
  labs(x = 'GBMI meta z (real)', y = 'z_synth at reported N',
       title = 'R10 Section 1 sanity — constructed vs real meta z',
       subtitle = 'chr22; identity line dashed red') +
  theme_half_open() + background_grid() +
  theme(strip.text = element_text(size = 10))
png(file.path(IMG, 'gbmi_r10_sanity.png'),
    res = 200, width = 2400, height = 1200, units = 'px')
print(p); dev.off()

# Gate: fail loud if any trait's r < threshold
fail <- res[r_zsynth_zmeta < THRESHOLD]
if (nrow(fail) > 0) {
  cat(sprintf('\nSANITY GATE FAILURE: %d traits below r=%.2f\n',
              nrow(fail), THRESHOLD))
  print(fail[, .(trait, r_zsynth_zmeta)])
  quit(status = 1)
}
cat(sprintf('\nAll %d traits pass r >= %.2f gate.\n', nrow(res), THRESHOLD))
