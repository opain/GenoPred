#!/usr/bin/env Rscript
# Round 7b Section 4 closing Task 1a — re-score OpenSNP EUR PGS restricted to
# the 6-way intersection of non-zero SNPs across all panels, then re-run the
# target R² comparison so that variant count is not a confounder.
#
# Strategy: use plink2 --score per panel per chr, applied to the shared
# OpenSNP target pgen, with --extract limited to the intersection SNP set.
# Sum per-sample per-chr sums into a genome-wide PGS. Then fit the same
# marginal r and paired-Δr bootstrap as m5_target_r2.R.

suppressPackageStartupMessages({
  library(data.table)
})

BASE   <- '/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred'
LABELS <- c('MELD_lambda0','EUR','EAS','AFR','CSA','AMR')
OUT_DIR <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
PHENO_PATH <- '/users/k1806347/oliverpainfel/Data/OpenSNP/processed/pheno/height.txt'
PLINK2 <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
B_BOOT <- 2000L
set.seed(2026L)

# ---------------------------------------------------------------------------
# Load intersection SNP set
snp_sets <- readRDS(file.path(OUT_DIR, 'meld_r7b_S4_close_snp_sets.rds'))
inter_all <- snp_sets$inter_all
cat(sprintf('Intersection SNPs across all 6 panels: %d\n', length(inter_all)))

# EUR keep from ancestry inference (any panel — same shared refdir)
eur_keep_path <- file.path(BASE, 'meld_test_r7bS4_ldpred2_MELD_lambda0',
                           'opensnp/ancestry/keep_files/model_based/EUR.keep')
if (!file.exists(eur_keep_path)) stop('EUR keep missing: ', eur_keep_path)

# ---------------------------------------------------------------------------
# Per panel: for each chr, run plink2 --score with --extract intersection.
# Concatenate across chr into a per-sample PGS.
tmp <- tempdir()
extract_path <- file.path(tmp, 'intersect.snps')
writeLines(inter_all, extract_path)

pgs_by_panel <- list()
for (L in LABELS) {
  cat(sprintf('scoring %s ... ', L))
  outdir <- file.path(BASE, sprintf('meld_test_r7bS4_ldpred2_%s', L))
  score  <- file.path(outdir, 'reference/pgs_score_files/ldpred2/yengo_all/ref-yengo_all.score.gz')
  # Merge chr-wise SSCORE sums
  per_chr_scores <- vector('list', 22L)
  for (chr in 1:22) {
    pgen <- file.path(outdir, sprintf('opensnp/geno/opensnp.ref.chr%d', chr))
    out  <- file.path(tmp, sprintf('score_%s_chr%d', L, chr))
    # plink2 --score: cols=+scoresums outputs the sum (not the mean)
    cmd <- c(PLINK2,
             '--pfile',   pgen,
             '--keep',    eur_keep_path,
             '--extract', extract_path,
             '--score',   score, '1', '2', '4', 'header',
                          'cols=+scoresums,-scoreavgs,-nallele,-dosagesum',
             '--out',     out)
    log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(log, 'status')) && attr(log, 'status') != 0L) {
      cat('\n'); cat(log, sep = '\n'); stop('plink2 --score failed for ', L, ' chr', chr)
    }
    d <- fread(paste0(out, '.sscore'))
    # Score sum column name after cols= tweak: "SCORE1_SUM"
    sc_col <- intersect(c('SCORE1_SUM','SCORE_SUM'), names(d))[1]
    if (is.na(sc_col)) sc_col <- grep('SUM$', names(d), value = TRUE)[1]
    per_chr_scores[[chr]] <- d[, .(FID = `#FID`, IID, chr_score = get(sc_col))]
  }
  # Sum across chr per (FID, IID)
  all_chr <- rbindlist(per_chr_scores)
  agg <- all_chr[, .(pgs_intersected = sum(chr_score)), by = .(FID, IID)]
  setnames(agg, 'pgs_intersected', L)
  pgs_by_panel[[L]] <- agg
  cat(sprintf(' %d samples\n', nrow(agg)))
}

# Merge all 6 panels on (FID, IID)
merged <- Reduce(function(a, b) merge(a, b, by = c('FID','IID')), pgs_by_panel)
cat(sprintf('\nMerged: %d samples with intersected PGS on all 6 panels\n', nrow(merged)))

# Merge with phenotype
pheno <- fread(PHENO_PATH)
d <- merge(pheno, merged, by = c('FID','IID'))
cat(sprintf('With height: %d individuals\n', nrow(d)))

# ---------------------------------------------------------------------------
# Bootstrap functions (same as m5_target_r2.R)
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

# ---------------------------------------------------------------------------
# Marginal r
marg <- rbindlist(lapply(LABELS, function(L) {
  st <- boot_r(d$height, d[[L]])
  data.table(panel = L,
             r     = st[['r']],
             r2    = st[['r']]^2,
             sd    = st[['sd']],
             lo95  = st[['lo95']], hi95 = st[['hi95']])
}))
marg[, r_txt  := sprintf('%.4f [%.4f, %.4f]', r, lo95, hi95)]
marg[, r2_txt := sprintf('%.4f', r2)]
setorder(marg, -r)
fwrite(marg, file.path(OUT_DIR, 'meld_r7b_S4_close_target_r2_intersected_marginal.csv'))

paired <- list()
for (L in setdiff(LABELS, 'MELD_lambda0')) {
  st <- boot_paired_diff_r(d$height, d$MELD_lambda0, d[[L]])
  paired[[length(paired) + 1L]] <- data.table(
    comparison = sprintf('MELD_lambda0 - %s', L),
    diff_r     = st[['diff']], sd = st[['sd']],
    lo95 = st[['lo95']], hi95 = st[['hi95']])
}
paired_dt <- rbindlist(paired)
paired_dt[, txt := sprintf('%+.4f [%+.4f, %+.4f]', diff_r, lo95, hi95)]
fwrite(paired_dt, file.path(OUT_DIR, 'meld_r7b_S4_close_target_r2_intersected_paired.csv'))

cat('\n=== INTERSECTED marginal r/r² (n=', nrow(d), ', ',
    length(inter_all), ' intersection SNPs) ===\n', sep = '')
print(marg[, .(panel, r_txt, r2_txt)])

cat('\n=== INTERSECTED paired Δr vs MELD-λ0 ===\n')
print(paired_dt[, .(comparison, txt)])

cat('\nDONE\n')
