#!/usr/bin/env Rscript
# Round 6, Section 1 helper.
# Per-population empirical dosage variance and ALT-allele frequency from the
# 1KG+HGDP reference. The Section 1 IVW meta target uses v_p (Section 1(v_p,i)),
# not the HWE 2f(1-f), so that residual within-superpopulation structure is
# absorbed into the weights.
#
# Usage:
#   Rscript build_ref_empirical_var.R <POP>       # POP ∈ {EUR, EAS, AFR, CSA, AMR, MID}
#
# Output:
#   misc/opensnp/ref_empirical/<POP>.var.rds   with columns
#     CHR, BP, SNP, REF, ALT, f_alt, v, n
#   (f_alt = mean(G_alt)/2, v = Var(G_alt) with n_called-1 denominator)

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: build_ref_empirical_var.R <POP>")
POP <- args[[1L]]

refdir     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
keep_file  <- file.path(refdir, 'keep_files', paste0(POP, '.keep'))
out_dir    <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
plink2     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'

stopifnot(file.exists(keep_file), file.exists(plink2))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tmp <- tempfile(pattern = paste0('empvar_', POP, '_'))
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

per_chr <- vector('list', 22L)
for (chr in 1:22) {
  pfile <- file.path(refdir, sprintf('ref.chr%d', chr))
  stopifnot(file.exists(paste0(pfile, '.pgen')))

  out_prefix <- file.path(tmp, sprintf('%s_chr%d', POP, chr))
  cmd <- c(plink2,
           '--pfile',        pfile,
           '--keep',         keep_file,
           '--geno-counts',
           '--out',          out_prefix)
  log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
  status <- attr(log, 'status')
  if (!is.null(status) && status != 0L) {
    cat(log, sep = '\n')
    stop(sprintf('plink2 geno-counts failed for POP=%s chr=%d', POP, chr))
  }

  gc  <- fread(paste0(out_prefix, '.gcount'))
  pv  <- fread(paste0(pfile, '.pvar'), select = c('#CHROM', 'POS', 'ID'))
  setnames(pv, c('#CHROM', 'POS', 'ID'), c('CHR', 'BP', 'ID'))
  setnames(gc, '#CHROM', 'CHR')

  # G = dosage of ALT (0/1/2). Counts (unphased diploid):
  #   HOM_REF_CT       = count of G=0
  #   HET_REF_ALT_CTS  = count of G=1
  #   TWO_ALT_GENO_CTS = count of G=2
  n <- gc$HOM_REF_CT + gc$HET_REF_ALT_CTS + gc$TWO_ALT_GENO_CTS
  sumG  <- gc$HET_REF_ALT_CTS + 2 * gc$TWO_ALT_GENO_CTS
  sumG2 <- gc$HET_REF_ALT_CTS + 4 * gc$TWO_ALT_GENO_CTS
  mean_G <- sumG / n
  # Population variance (n_called denominator) then Bessel-correct to sample var
  var_G_pop  <- sumG2 / n - mean_G^2
  # sample variance = pop_variance * n / (n-1)
  var_G_samp <- var_G_pop * n / pmax(n - 1L, 1L)
  # Handle SNPs with n < 2 (no variance): set NA and drop downstream.
  var_G_samp[n < 2L] <- NA_real_

  out <- data.table(
    CHR   = gc$CHR,
    BP    = pv$BP[match(gc$ID, pv$ID)],
    SNP   = gc$ID,
    REF   = gc$REF,
    ALT   = gc$ALT,
    f_alt = mean_G / 2,
    v     = var_G_samp,
    n     = n
  )
  per_chr[[chr]] <- out
  cat(sprintf('POP=%s chr=%d: %d SNPs, %d n<2\n', POP, chr, nrow(out), sum(n < 2L)))
}

all <- rbindlist(per_chr)
saveRDS(all, file.path(out_dir, paste0(POP, '.var.rds')))
cat(sprintf('wrote %s (%d SNPs)\n',
            file.path(out_dir, paste0(POP, '.var.rds')), nrow(all)))
