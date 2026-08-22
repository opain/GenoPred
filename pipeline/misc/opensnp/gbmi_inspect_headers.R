#!/usr/bin/env Rscript
# Inspect GBMI file headers to derive the column mapping for harmonisation.
# Reads a small sample from ONE meta file and ONE arm file and reports:
#   - column names
#   - dtypes (from data.table auto-typing)
#   - genome-build indicator: chromosome codes + first-position sample
#   - presence of AF (allele frequency), per-SNP N, effect (beta/OR), SE columns
#   - presence of per-biobank or per-ancestry columns
# Output: gbmi_r8_header_report.txt + prints to stdout.

suppressPackageStartupMessages(library(data.table))

DEST     <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI'
META_FN  <- 'Asthma_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'
ARM_FN   <- 'Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'
OUT      <- file.path(DEST, 'gbmi_r8_header_report.txt')

sink(OUT, split = TRUE)
on.exit(sink())

inspect <- function(fp, label) {
  cat(sprintf('\n\n===== %s =====\n', label))
  cat(sprintf('file: %s\n', fp))
  cat(sprintf('size: %.2f GB\n', file.info(fp)$size / 1024^3))

  # Read header + first 20 000 rows for column typing and quick sanity
  d <- fread(fp, nrows = 20000, showProgress = FALSE)
  cat(sprintf('rows read: %d\n', nrow(d)))
  cat(sprintf('n cols: %d\n', ncol(d)))
  cat('\ncolumn names + types:\n')
  print(data.table(col = names(d), type = sapply(d, function(x) class(x)[1])))

  cat('\nfirst 5 rows:\n')
  print(head(d, 5))

  cat('\nchromosome codes (distribution):\n')
  chr_col <- intersect(names(d), c('CHR','#CHR','chr','chrom','CHROM','#CHROM'))
  if (length(chr_col) > 0) {
    cat(sprintf('  chr column: %s\n', chr_col[1]))
    print(table(d[[chr_col[1]]]))
  } else cat('  no obvious chr column\n')

  cat('\nposition column diagnostics (min/max for build detection):\n')
  pos_col <- intersect(names(d), c('POS','BP','bp','pos','position','base_pair_location'))
  if (length(pos_col) > 0) {
    cat(sprintf('  pos column: %s\n', pos_col[1]))
    cat(sprintf('  min/max: %d / %d\n', min(d[[pos_col[1]]]), max(d[[pos_col[1]]])))
  } else cat('  no obvious pos column\n')

  cat('\nlikely-effect / SE / P columns (name grep):\n')
  eff_cols <- names(d)[grepl('beta|effect|OR|logOR', names(d), ignore.case = TRUE)]
  se_cols  <- names(d)[grepl('^se|_se|standard.?error', names(d), ignore.case = TRUE)]
  p_cols   <- names(d)[grepl('^p[_-]?val|^p$|pvalue', names(d), ignore.case = TRUE)]
  af_cols  <- names(d)[grepl('af|freq|maf', names(d), ignore.case = TRUE)]
  n_cols   <- names(d)[grepl('^n$|^n[_-]|sample.?size|neff', names(d), ignore.case = TRUE)]
  cat(sprintf('  effect: %s\n', paste(eff_cols, collapse = ', ')))
  cat(sprintf('  SE: %s\n', paste(se_cols, collapse = ', ')))
  cat(sprintf('  P: %s\n', paste(p_cols, collapse = ', ')))
  cat(sprintf('  AF: %s\n', paste(af_cols, collapse = ', ')))
  cat(sprintf('  N: %s\n', paste(n_cols, collapse = ', ')))

  cat('\nper-biobank / per-ancestry columns (name grep):\n')
  bb_cols <- names(d)[grepl('BioMe|BioVU|CCPM|DECODE|ESTBB|FinnGen|HUNT|MGB|MGI|UKBB|biobank', names(d), ignore.case = TRUE)]
  anc_cols <- names(d)[grepl('_eur|_afr|_amr|_eas|_sas|_mid|ancestry', names(d), ignore.case = TRUE)]
  cat(sprintf('  biobank: %s\n', paste(bb_cols, collapse = ', ')))
  cat(sprintf('  ancestry: %s\n', paste(anc_cols, collapse = ', ')))
  invisible(d)
}

inspect(file.path(DEST, META_FN), 'META (Asthma Bothsex all-ancestry)')
inspect(file.path(DEST, ARM_FN),  'ARM (Asthma Bothsex eur)')

cat('\n\n===== notes =====\n')
cat('GBMI files as documented in Zhou et al. Cell Genomics 2022 supplementary:\n')
cat('  - Build GRCh38 (b38) — needs confirmation from actual position ranges.\n')
cat('  - Columns typically: CHR POS REF ALT rsid all_meta_N all_meta_AF all_inv_var_meta_beta\n')
cat('    all_inv_var_meta_sebeta all_inv_var_meta_p all_inv_var_het_p ...plus per-biobank blocks.\n')
cat('  Meta files also carry per-ancestry columns for each of the ancestries listed in the manifest.\n')
