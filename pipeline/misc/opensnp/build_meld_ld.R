#!/usr/bin/env Rscript
# Round 6, Section 2. Build per-LDetect-block empirical LD matrices and the
# analytic ground-truth ingredients (R_p, v_p, f_p) needed to reconstruct the
# meta target R* = D*^{-1/2} Σ* D*^{-1/2} and the mega target R_pool at any
# mixture weight w without recomputing correlations.
#
# For each block, saves an RDS with:
#   SNP, A1, A2, BP, cM, f_EUR, f_AFR, v_EUR, v_AFR, R_EUR, R_AFR
# where R_p is the dense within-block correlation matrix computed with the
# same 3 cM bigsnpr window as GenoPred's LDpred2 rules use.
#
# Runs the block-wise invariants inline; a fuller cross-check vs empirical
# pooled LD lives in tests/test_build_meld_ld.R (invariant 4).
#
# Usage:
#   Rscript build_meld_ld.R [<chr_from> [<chr_to>]]        # defaults: 1..22

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# Config
refdir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
plink2      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
map_dir     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_map/genetic_map_b37'
block_dir   <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ld_blocks/EUR'
emp_var_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
out_dir     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
pops        <- c('EUR', 'AFR')
window_cm   <- 3                      # matches GenoPred LDpred2 (size = 3/1000 M = 3 cM) and quickprs LDAK --window-cm 3
# Use ncores = 1 for snp_cor: bigsnpr complains about nested parallelism if BLAS
# is already multithreaded, and per-block problem sizes (few hundred SNPs × ~700
# samples) don't benefit meaningfully from >1 thread.
ncores      <- 1

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args      <- commandArgs(trailingOnly = TRUE)
chr_from  <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
chr_to    <- if (length(args) >= 2) as.integer(args[[2]]) else 22L

# ---------------------------------------------------------------------------
# Helpers
interp_cM <- function(bp, map_bp, map_cM) {
  # Piecewise-linear interpolation with extrapolation to endpoints via nearest cM.
  approx(map_bp, map_cM, xout = bp, rule = 2L)$y
}

read_block_bed <- function(chr) {
  bf <- fread(file.path(block_dir, sprintf('fourier_ls-chr%d.bed', chr)))
  setnames(bf, tolower(trimws(names(bf))))
  bf$start <- as.integer(bf$start)
  bf$stop  <- as.integer(bf$stop)
  bf[order(start)]
}

read_pop_ids <- function(pop) {
  readLines(file.path(refdir, 'keep_files', paste0(pop, '.keep')))
}

read_genmap <- function(chr) {
  fread(file.path(map_dir, sprintf('genetic_map_chr%d_combined_b37.txt', chr)))
}

# ---------------------------------------------------------------------------
# Load empirical variance/frequency once (used for f_p, v_p; also gives us
# the canonical SNP set and allele orientation that the sumstats agree with)
emp <- list()
for (p in pops) emp[[p]] <- readRDS(file.path(emp_var_dir, paste0(p, '.var.rds')))
stopifnot(all(emp$EUR$SNP == emp$AFR$SNP))

# Per-population sample-id → row index (built once per chr from bigSNP fam)
pop_ids <- lapply(setNames(pops, pops), read_pop_ids)

# ---------------------------------------------------------------------------
# Per-chromosome outer loop
for (chr in chr_from:chr_to) {
  cat(sprintf('\n=== chromosome %d ===\n', chr))
  chr_out_dir <- file.path(out_dir, sprintf('chr%d', chr))
  dir.create(chr_out_dir, showWarnings = FALSE, recursive = TRUE)

  # (a) Convert pgen → bed once (temp dir; auto-cleanup on exit)
  tmp <- tempfile(pattern = sprintf('mldld_chr%d_', chr))
  dir.create(tmp)
  bed_prefix <- file.path(tmp, sprintf('chr%d', chr))
  cmd <- c(plink2,
           '--pfile',      file.path(refdir, sprintf('ref.chr%d', chr)),
           '--make-bed',
           '--out',        bed_prefix)
  log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
  status <- attr(log, 'status')
  if (!is.null(status) && status != 0L) {
    cat(log, sep = '\n')
    unlink(tmp, recursive = TRUE)
    stop(sprintf('plink2 --make-bed failed for chr=%d', chr))
  }

  # (b) Attach bigSNP
  bk <- snp_readBed2(paste0(bed_prefix, '.bed'),
                     backingfile = tempfile(pattern = sprintf('bk_chr%d_', chr)))
  obj <- snp_attach(bk)
  G   <- obj$genotypes
  fam <- obj$fam
  bim <- obj$map
  setDT(bim); setDT(fam)

  # (c) Interpolate cM per SNP
  gm <- read_genmap(chr)
  bim[, cM := interp_cM(bim$physical.pos, gm$position, gm$`Genetic_Map(cM)`)]

  # (d) Sample indices per pop
  fam[, row := .I]
  ind_row_by_pop <- lapply(setNames(pops, pops), function(p) {
    idx <- match(pop_ids[[p]], fam$sample.ID)
    if (any(is.na(idx))) stop(sprintf('unmatched IIDs for pop=%s (%d missing)',
                                      p, sum(is.na(idx))))
    idx
  })

  # (e) Read block bed
  bf <- read_block_bed(chr)
  cat(sprintf('  %d blocks; %d SNPs on chr%d\n', nrow(bf), nrow(bim), chr))

  # (f) Per-block loop
  for (bi in seq_len(nrow(bf))) {
    b_start <- bf$start[bi]
    b_stop  <- bf$stop[bi]
    # LDetect blocks are half-open [start, stop); include SNPs at start, exclude stop
    idx_col <- bim[physical.pos >= b_start & physical.pos < b_stop, which = TRUE]
    if (length(idx_col) < 2L) next        # can't correlate <2 SNPs

    # Match empirical variance/AF onto these SNPs, in the same allele orientation.
    # emp uses REF/ALT from the pvar; bim contains allele1/allele2 which are A1/A2
    # in the .bim (plink1 convention). We need to align so that v/f match the
    # dosage direction used by snp_cor. bigsnpr snp_cor uses G = dosage of A1
    # (allele1) in the bim — which after plink2 --make-bed defaults to ALT.
    # Confirm:
    block_snps <- bim$marker.ID[idx_col]
    m <- match(block_snps, emp$EUR$SNP)
    if (any(is.na(m))) {
      # SNPs in bed but not in emp — skip those (rare edge case)
      keep_local <- !is.na(m)
      idx_col   <- idx_col[keep_local]
      block_snps <- bim$marker.ID[idx_col]
      m <- match(block_snps, emp$EUR$SNP)
    }
    if (length(idx_col) < 2L) next

    # Allele orientation check: in the bim from plink2 --make-bed, allele1
    # ("A1") is REF and allele2 ("A2") is ALT (plink2 default swap of the plink1
    # convention). Confirm and record; if this ever surprises us we should
    # fail loudly rather than silently mis-orient.
    a1_bim <- bim$allele1[idx_col]
    a2_bim <- bim$allele2[idx_col]
    ref_e  <- emp$EUR$REF[m]
    alt_e  <- emp$EUR$ALT[m]

    same <- (a1_bim == ref_e & a2_bim == alt_e)
    flip <- (a1_bim == alt_e & a2_bim == ref_e)
    if (!all(same | flip)) {
      stop(sprintf('chr%d block%d: bim/emp allele mismatch on %d SNPs',
                   chr, bi, sum(!(same | flip))))
    }

    # We want R and v aligned to the ALT-coded canonical direction. v_p is
    # invariant to allele swap, so nothing to do there. The DIAGONAL of R_p is
    # invariant too (=1). But pairwise off-diagonal correlations Cor(G_i, G_j)
    # DO flip sign if i and j are coded on opposite strands: swap A1 on one SNP
    # → G → 2−G, and Cor(2−G_i, G_j) = −Cor(G_i, G_j). So if bim mixes REF-coded
    # and ALT-coded SNPs, we must sign-flip R[i,j] where sign(i) ≠ sign(j).
    #    same → R_p was computed on G_p = dose of REF (bim A1 = REF, our target
    #           is ALT-coded → sign = −1)
    #    flip → R_p was computed on G_p = dose of ALT (bim A1 = ALT → sign = +1)
    # With plink2 --make-bed the default is A1=ALT so flip covers all SNPs; but
    # ref-major SNPs or edge cases can produce mixed orientation, so make the
    # correction unconditional.
    sign_canon <- ifelse(flip, +1, -1)     # +1 if bim A1 == ALT, −1 if bim A1 == REF
    # f_p is stored as f_alt (ALT-coded) so already matches canonical — no flip.

    # Genetic positions in Morgans (snp_cor size convention in GenoPred is 3/1000
    # in Morgans; equivalently 3 cM). Reproduce that idiom exactly.
    pos_M <- bim$cM[idx_col] / 100
    size_M <- window_cm / 100

    R_list <- list(); v_list <- list(); f_list <- list()
    for (p in pops) {
      R_p_sparse <- snp_cor(G,
                            ind.row   = ind_row_by_pop[[p]],
                            ind.col   = idx_col,
                            size      = size_M,
                            infos.pos = pos_M,
                            ncores    = ncores)
      R_p <- as.matrix(R_p_sparse)
      # snp_cor may return NaN diagonals for zero-variance SNPs (monomorphic).
      diag(R_p)[is.na(diag(R_p))] <- 1
      R_p[is.na(R_p)] <- 0
      # Sign-flip to canonical ALT-coded orientation (see comment above).
      R_p <- R_p * outer(sign_canon, sign_canon)
      R_list[[p]] <- R_p

      v_list[[p]] <- emp[[p]]$v[m]
      f_list[[p]] <- emp[[p]]$f_alt[m]
    }

    # Invariant 1 spot-check (per-block): diag(R_p) == 1 to tolerance
    for (p in pops) {
      dd <- diag(R_list[[p]])
      if (max(abs(dd - 1)) > 1e-8 & any(v_list[[p]] > 0)) {
        stop(sprintf('chr%d block%d %s: max|diag-1| = %.3e (v>0 SNPs present)',
                     chr, bi, p, max(abs(dd - 1))))
      }
    }

    out_bl <- list(
      chr        = chr,
      block_id   = bi,
      block_bp   = c(b_start, b_stop),
      SNP        = block_snps,
      BP         = bim$physical.pos[idx_col],
      cM         = bim$cM[idx_col],
      A1_canon   = alt_e,                # canonical A1 = reference ALT
      A2_canon   = ref_e,
      f_EUR      = f_list$EUR,
      f_AFR      = f_list$AFR,
      v_EUR      = v_list$EUR,
      v_AFR      = v_list$AFR,
      R_EUR      = R_list$EUR,
      R_AFR      = R_list$AFR
    )
    saveRDS(out_bl,
            file.path(chr_out_dir, sprintf('block_%04d.rds', bi)),
            compress = 'gzip')
  }
  cat(sprintf('  wrote %d block files under %s\n',
              length(list.files(chr_out_dir, pattern = '^block_.*rds$')),
              chr_out_dir))

  # Cleanup bigsnpr backing files for this chr
  bk_files <- c(bk, sub('\\.bk$', '.rds', bk))
  suppressWarnings(file.remove(bk_files[file.exists(bk_files)]))
  unlink(tmp, recursive = TRUE)
}

cat('\nDONE\n')
