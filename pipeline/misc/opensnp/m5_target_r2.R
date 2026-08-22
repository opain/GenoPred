#!/usr/bin/env Rscript
# Round 7b Section 4 — target R² per panel.
# For each of the 6 LD panels, fit lm(scale(height) ~ scale(PGS)) on OpenSNP
# EUR-inferred samples, report r, r², and bootstrap 95% CIs (resample
# individuals with replacement). Also compute a paired bootstrap of (r_A − r_B)
# for MELD-λ0 vs each single-pop, to control for shared sampling noise.

suppressPackageStartupMessages(library(data.table))

BASE   <- '/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred'
LABELS <- c('MELD_lambda0','EUR','EAS','AFR','CSA','AMR','equal')
OUT_DIR <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
PHENO_PATH <- '/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno/height.txt'
B_BOOT <- 2000L
set.seed(2026L)

pheno <- fread(PHENO_PATH)
cat(sprintf('Phenotype rows: %d\n', nrow(pheno)))

# Load PGS from each panel's OpenSNP-EUR target scoring
panels <- list()
for (L in LABELS) {
  fp <- file.path(BASE, sprintf('meld_test_r7bS4_ldpred2_%s', L),
                  'opensnp/pgs/EUR/ldpred2/yengo_all/opensnp-yengo_all-EUR.raw.profiles')
  d <- fread(fp)
  sc <- setdiff(names(d), c('FID','IID'))[1L]
  d <- d[, .(FID, IID, pgs = get(sc))]
  setnames(d, 'pgs', L)
  panels[[L]] <- d
}
merged <- Reduce(function(a, b) merge(a, b, by = c('FID','IID')), panels)
d <- merge(pheno, merged, by = c('FID','IID'))
cat(sprintf('OpenSNP EUR with height + PGS on all 6 panels: %d individuals\n', nrow(d)))

# Marginal r per panel with bootstrap CI
boot_r <- function(y, x, B = B_BOOT) {
  ok <- is.finite(y) & is.finite(x)
  n <- sum(ok); yo <- y[ok]; xo <- x[ok]
  r0 <- suppressWarnings(cor(yo, xo))
  boots <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(yo[i], xo[i]))
  })
  c(r = r0, sd = sd(boots), lo95 = as.numeric(quantile(boots, 0.025)),
    hi95 = as.numeric(quantile(boots, 0.975)))
}

marg <- rbindlist(lapply(LABELS, function(L) {
  st <- boot_r(d$height, d[[L]])
  data.table(panel = L,
             r     = st[['r']],
             r2    = st[['r']]^2,
             sd    = st[['sd']],
             lo95  = st[['lo95']], hi95 = st[['hi95']])
}))
marg[, r_txt   := sprintf('%.4f [%.4f, %.4f]', r, lo95, hi95)]
marg[, r2_txt  := sprintf('%.4f', r2)]
setorder(marg, -r)
fwrite(marg, file.path(OUT_DIR, 'meld_r7b_S4_target_r2_marginal.csv'))

# Paired difference (r_A − r_B) with paired bootstrap — same resample of individuals
boot_paired_diff_r <- function(y, xa, xb, B = B_BOOT) {
  ok <- is.finite(y) & is.finite(xa) & is.finite(xb)
  n <- sum(ok); yo <- y[ok]; xao <- xa[ok]; xbo <- xb[ok]
  d0 <- suppressWarnings(cor(yo, xao) - cor(yo, xbo))
  boots <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(yo[i], xao[i]) - cor(yo[i], xbo[i]))
  })
  c(diff = d0, sd = sd(boots),
    lo95 = as.numeric(quantile(boots, 0.025)),
    hi95 = as.numeric(quantile(boots, 0.975)))
}

paired <- list()
for (L in setdiff(LABELS, 'MELD_lambda0')) {
  st <- boot_paired_diff_r(d$height, d$MELD_lambda0, d[[L]])
  paired[[length(paired) + 1L]] <- data.table(
    comparison = sprintf('MELD_lambda0 − %s', L),
    diff_r     = st[['diff']], sd = st[['sd']],
    lo95 = st[['lo95']], hi95 = st[['hi95']])
}
# R7b S4 closing T2: also report equal − EUR to isolate composition from shrinkage.
if (all(c('equal','EUR') %in% LABELS)) {
  st <- boot_paired_diff_r(d$height, d$equal, d$EUR)
  paired[[length(paired) + 1L]] <- data.table(
    comparison = 'equal − EUR',
    diff_r     = st[['diff']], sd = st[['sd']],
    lo95 = st[['lo95']], hi95 = st[['hi95']])
}
paired_dt <- rbindlist(paired)
paired_dt[, txt := sprintf('%+.4f [%+.4f, %+.4f]', diff_r, lo95, hi95)]
fwrite(paired_dt, file.path(OUT_DIR, 'meld_r7b_S4_target_r2_paired_diff.csv'))

cat('\n=== Marginal r and r² (OpenSNP EUR height ~ PGS, n=', nrow(d), ') ===\n', sep = '')
print(marg[, .(panel, r_txt, r2_txt)])

cat('\n=== Paired Δr vs MELD-λ0 (95% bootstrap CI) ===\n')
print(paired_dt[, .(comparison, txt)])

cat('\nDONE\n')
