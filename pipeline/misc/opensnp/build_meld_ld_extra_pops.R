#!/usr/bin/env Rscript
# Round 6 Section 3 helper: extend the per-block LD RDS files produced by
# build_meld_ld.R with additional single-population correlation matrices
# (default EAS, CSA, AMR) so that eval_meld_ld_target_free.R can compare
# every single-pop candidate against the analytic meta target R* without
# re-attaching pgen files at eval time.
#
# For each chromosome + block, computes snp_cor on the block SNPs restricted
# to the requested pop's keep_file, sign-aligns to the canonical ALT allele,
# and appends fields R_<POP> and v_<POP> to the existing block RDS.
#
# Idempotent: rewrites each block file only if the requested pop's field is
# missing. Skips blocks that already have all requested pops.
#
# Usage:
#   Rscript build_meld_ld_extra_pops.R [<POP1[,POP2,...]>] [<chr_from> [<chr_to>]]
#     default POPs: EAS,CSA,AMR

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
refdir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
plink2      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
emp_var_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
meld_ld_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
window_cm   <- 3

args     <- commandArgs(trailingOnly = TRUE)
pops     <- if (length(args) >= 1) strsplit(args[[1]], ',')[[1]] else c('EAS','CSA','AMR')
chr_from <- if (length(args) >= 2) as.integer(args[[2]]) else 1L
chr_to   <- if (length(args) >= 3) as.integer(args[[3]]) else 22L
stopifnot(all(pops %in% c('EUR','EAS','AFR','CSA','AMR','MID')))

pop_ids  <- lapply(setNames(pops, pops), function(p)
  readLines(file.path(refdir, 'keep_files', paste0(p, '.keep'))))
emp      <- lapply(setNames(pops, pops), function(p)
  readRDS(file.path(emp_var_dir, paste0(p, '.var.rds'))))

for (chr in chr_from:chr_to) {
  chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', chr))
  block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
  if (length(block_files) == 0L) {
    cat(sprintf('chr%d: no block files, skipping\n', chr)); next
  }

  # Check which blocks still need work
  need <- FALSE
  for (fn in block_files) {
    bl <- readRDS(fn)
    if (any(!(paste0('R_', pops) %in% names(bl)))) { need <- TRUE; break }
  }
  if (!need) { cat(sprintf('chr%d: all blocks already have %s — skip\n',
                           chr, paste(pops, collapse=','))); next }

  cat(sprintf('\n=== chr%d: extending %d block files with %s ===\n',
              chr, length(block_files), paste(pops, collapse=',')))

  # Convert pgen → bed once
  tmp <- tempfile(pattern = sprintf('mldext_chr%d_', chr))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  bed_prefix <- file.path(tmp, sprintf('chr%d', chr))
  cmd <- c(plink2, '--pfile', file.path(refdir, sprintf('ref.chr%d', chr)),
           '--make-bed', '--out', bed_prefix)
  log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(log, 'status'))) { cat(log, sep='\n'); stop('plink2 make-bed failed') }

  bk  <- snp_readBed2(paste0(bed_prefix, '.bed'),
                      backingfile = tempfile(pattern = sprintf('bkext_chr%d_', chr)))
  obj <- snp_attach(bk)
  G   <- obj$genotypes
  fam <- as.data.table(obj$fam)
  bim <- as.data.table(obj$map)
  fam[, row := .I]

  ind_row_by_pop <- lapply(setNames(pops, pops), function(p) {
    idx <- match(pop_ids[[p]], fam$sample.ID)
    if (any(is.na(idx))) stop(sprintf('unmatched IIDs for pop=%s (%d missing)',
                                      p, sum(is.na(idx))))
    idx
  })

  written <- 0L
  for (fn in block_files) {
    bl <- readRDS(fn)
    missing_pops <- pops[!paste0('R_', pops) %in% names(bl)]
    if (length(missing_pops) == 0L) next

    # Find SNP indices in the bim
    ord <- match(bl$SNP, bim$marker.ID)
    ok  <- !is.na(ord)
    if (sum(ok) < 2L) next

    a1_bim <- bim$allele1[ord[ok]]
    a2_bim <- bim$allele2[ord[ok]]
    alt_e  <- bl$A1_canon[ok]
    ref_e  <- bl$A2_canon[ok]
    same   <- (a1_bim == ref_e & a2_bim == alt_e)
    flip   <- (a1_bim == alt_e & a2_bim == ref_e)
    if (!all(same | flip)) stop(sprintf('chr%d block%d: bim/canon allele mismatch',
                                        bl$chr, bl$block_id))
    sign_canon <- ifelse(flip, +1, -1)
    pos_M <- bl$cM[ok] / 100

    for (p in missing_pops) {
      R_sp <- snp_cor(G, ind.row = ind_row_by_pop[[p]],
                      ind.col = ord[ok],
                      size = window_cm / 100,
                      infos.pos = pos_M, ncores = 1)
      R <- as.matrix(R_sp)
      diag(R)[is.na(diag(R))] <- 1
      R[is.na(R)] <- 0
      R <- R * outer(sign_canon, sign_canon)
      # Pad R back to full block SNP set (fill zeros for dropped rows/cols)
      if (!all(ok)) {
        R_full <- matrix(0, length(bl$SNP), length(bl$SNP))
        R_full[ok, ok] <- R
        diag(R_full) <- 1                 # keep diag = 1 for the padded rows too
        R <- R_full
      }
      m_emp <- match(bl$SNP, emp[[p]]$SNP)
      v_p <- emp[[p]]$v[m_emp];  v_p[is.na(v_p)] <- 0
      f_p <- emp[[p]]$f_alt[m_emp]

      bl[[paste0('R_', p)]] <- R
      bl[[paste0('v_', p)]] <- v_p
      bl[[paste0('f_', p)]] <- f_p
    }

    saveRDS(bl, fn, compress = 'gzip')
    written <- written + 1L
  }
  cat(sprintf('  extended %d block files on chr%d\n', written, chr))

  # cleanup backing files for this chr
  bk_files <- c(bk, sub('\\.bk$', '.rds', bk))
  suppressWarnings(file.remove(bk_files[file.exists(bk_files)]))
}

cat('\nDONE\n')
