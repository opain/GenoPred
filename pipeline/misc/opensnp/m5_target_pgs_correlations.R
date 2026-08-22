#!/usr/bin/env Rscript
# Round 7b Section 4 M5 (target-PGS version):
# For each of the 6 candidate LD panels, read the OpenSNP EUR target-PGS
# profile (LDpred2-auto beta), align samples across panels, and compute
# pairwise Pearson correlations of PGS values across panels.
#
# Outputs (pipeline/misc/opensnp/):
#   meld_r7b_S4_target_pgs_pairwise.csv       pairwise r + bootstrap CIs
#   meld_r7b_S4_target_pgs_matrix.csv         6×6 matrix of pairwise r
#   meld_r7b_S4_target_pgs_summary.csv        each panel vs MELD-λ0

suppressPackageStartupMessages(library(data.table))

BASE  <- '/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred'
LABELS <- c('MELD_lambda0','EUR','EAS','AFR','CSA','AMR','equal')
OUT_DIR <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
B_BOOT <- 2000L
set.seed(2026L)

# Load each panel's EUR profile
panels <- list()
for (L in LABELS) {
  fp <- file.path(BASE, sprintf('meld_test_r7bS4_ldpred2_%s', L),
                  'opensnp/pgs/EUR/ldpred2/yengo_all/opensnp-yengo_all-EUR.raw.profiles')
  if (!file.exists(fp)) stop('missing: ', fp)
  d <- fread(fp)
  # Score column name — typically yengo_all_beta_auto
  sc <- setdiff(names(d), c('FID','IID'))[1L]
  d <- d[, .(FID, IID, pgs = get(sc))]
  setnames(d, 'pgs', L)
  panels[[L]] <- d
  cat(sprintf('%s: %d samples (score col: %s)\n', L, nrow(d), sc))
}

# Align on (FID, IID)
merged <- Reduce(function(a, b) merge(a, b, by = c('FID','IID')), panels)
cat(sprintf('\nmerged: %d samples with PGS on all %d panels\n', nrow(merged), length(LABELS)))

M <- as.matrix(merged[, ..LABELS])

# Pairwise Pearson r + bootstrap CI (SNP-level bootstrap: resample individuals)
boot_pair <- function(x, y, B = B_BOOT) {
  ok <- is.finite(x) & is.finite(y)
  r0 <- suppressWarnings(cor(x[ok], y[ok]))
  n <- sum(ok); xo <- x[ok]; yo <- y[ok]
  boots <- replicate(B, {
    i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(xo[i], yo[i]))
  })
  c(r = r0, sd = sd(boots, na.rm = TRUE),
    lo95 = as.numeric(quantile(boots, 0.025, na.rm = TRUE)),
    hi95 = as.numeric(quantile(boots, 0.975, na.rm = TRUE)))
}

pair_rows <- list()
for (i in 1:(length(LABELS) - 1L)) {
  for (j in (i + 1L):length(LABELS)) {
    a <- LABELS[i]; b <- LABELS[j]
    st <- boot_pair(M[, a], M[, b])
    pair_rows[[length(pair_rows) + 1L]] <- data.table(
      panel_a = a, panel_b = b,
      r = st[['r']], sd = st[['sd']],
      lo95 = st[['lo95']], hi95 = st[['hi95']]
    )
  }
}
pair_dt <- rbindlist(pair_rows)
pair_dt[, txt := sprintf('%.4f [%.4f, %.4f]', r, lo95, hi95)]
fwrite(pair_dt, file.path(OUT_DIR, 'meld_r7b_S4_target_pgs_pairwise.csv'))

# Matrix form
mat <- matrix(1, length(LABELS), length(LABELS), dimnames = list(LABELS, LABELS))
for (i in seq_len(nrow(pair_dt))) {
  a <- pair_dt$panel_a[i]; b <- pair_dt$panel_b[i]
  mat[a, b] <- pair_dt$r[i]; mat[b, a] <- pair_dt$r[i]
}
fwrite(data.table(panel = LABELS, mat), file.path(OUT_DIR, 'meld_r7b_S4_target_pgs_matrix.csv'))

# Summary: each panel vs MELD-λ0
summary_dt <- pair_dt[panel_a == 'MELD_lambda0' | panel_b == 'MELD_lambda0']
summary_dt[, other := ifelse(panel_a == 'MELD_lambda0', panel_b, panel_a)]
summary_dt <- summary_dt[, .(other, r, sd, lo95, hi95, txt)]
setorder(summary_dt, -r)
fwrite(summary_dt, file.path(OUT_DIR, 'meld_r7b_S4_target_pgs_summary.csv'))

cat('\n=== Pairwise Pearson r of target PGS across 6 panels (bootstrap 95% CI) ===\n')
print(pair_dt[, .(panel_a, panel_b, txt)])
cat('\n=== 6×6 correlation matrix ===\n')
print(round(mat, 4))
cat('\n=== Each panel vs MELD-λ0 (sorted by r) ===\n')
print(summary_dt[, .(other, txt)])
cat('\nDONE\n')
