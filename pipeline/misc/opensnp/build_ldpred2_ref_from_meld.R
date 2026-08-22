#!/usr/bin/env Rscript
# Round 7b Section 4 — build an LDpred2-consumable reference directory from the
# analytic MELD-λ0 R* target (or a single-population R). Writes:
#   <out_dir>/LD_with_blocks_chr{1..22}.rds   as Matrix::dsCMatrix
#   <out_dir>/map.rds                         data.frame with chr, pos, a0, a1, af, ld, group_id
#
# Per LDetect block:
#   1. Reconstruct R for the chosen target (single-pop R_<pop> or MELD-λ0 R*)
#   2. Apply uniform eigen-clamp PSD repair (from Round 6b) so LDpred2's Gibbs
#      sampler stays stable — the M2 numerical-conditioning bug is not
#      documented to affect LDpred2 but the repair costs nothing.
#   3. Stitch per-block dense sub-blocks into a chr-wide block-diagonal
#      sparse dsCMatrix.
#
# Usage:
#   Rscript build_ldpred2_ref_from_meld.R <out_dir> <label>
# where <label> ∈ {'MELD_lambda0', 'EUR', 'EAS', 'AFR', 'CSA', 'AMR'}.
# For MELD_lambda0 the weights are AF-projection weights from Section 2
# (EUR 0.7800, EAS 0.1274, AFR 0.0528, CSA 0.0193, AMR 0.0205), hard-coded here
# to match the Section-2 headline exactly.

suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop('usage: build_ldpred2_ref_from_meld.R <out_dir> <label>')
OUT_DIR <- args[[1L]]
LABEL   <- args[[2L]]
stopifnot(LABEL %in% c('MELD_lambda0','EUR','EAS','AFR','CSA','AMR','equal'))

# ---------------------------------------------------------------------------
meld_ld_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
map_dir     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_map/genetic_map_b37'
pops        <- c('EUR','EAS','AFR','CSA','AMR')
# Composition weights for MELD-flavour labels. MELD_lambda0 uses AF-projection
# weights from Section 2 (78/13/5/2/2). 'equal' uses 0.2 per population; it is
# the isolate-shrinkage-from-composition control panel for R7b S4 closing T2.
w_afproj_meld  <- c(EUR = 0.7800, EAS = 0.1274, AFR = 0.0528, CSA = 0.0193, AMR = 0.0205)
w_afproj_meld  <- w_afproj_meld / sum(w_afproj_meld)
w_afproj_equal <- setNames(rep(1/length(pops), length(pops)), pops)

w_afproj <- if (LABEL == 'equal') w_afproj_equal else w_afproj_meld

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf('label=%s\nout_dir=%s\n', LABEL, OUT_DIR))
if (LABEL %in% c('MELD_lambda0','equal')) { cat('weights:\n'); print(round(w_afproj, 4)) }

# ---------------------------------------------------------------------------
# Helpers
psd_repair <- function(R) {
  ev <- eigen(R, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  R2 <- ev$vectors %*% (lam * t(ev$vectors))
  d <- sqrt(pmax(diag(R2), 0)); d[d == 0] <- 1
  R2 <- R2 / outer(d, d); diag(R2) <- 1
  R2
}

reconstruct_block <- function(bl, label, w = w_afproj) {
  if (!(label %in% c('MELD_lambda0','equal'))) {
    R <- bl[[paste0('R_', label)]]
    R[!is.finite(R)] <- 0
    return(R)
  }
  # MELD-flavour: P-pop meta target R*
  Sigma_star <- matrix(0, length(bl$SNP), length(bl$SNP))
  V_star     <- rep(0, length(bl$SNP))
  for (p in pops) {
    R_p <- bl[[paste0('R_', p)]]
    v_p <- bl[[paste0('v_', p)]]
    Cov_p <- R_p * outer(sqrt(v_p), sqrt(v_p))
    Sigma_star <- Sigma_star + w[[p]] * Cov_p
    V_star     <- V_star     + w[[p]] * v_p
  }
  inv_sd <- 1 / sqrt(V_star)
  R <- Sigma_star * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

# Per-SNP composition-weighted AF and per-SNP LD score (sum of r² within block)
build_map_row <- function(bl, R_rep, label, w = w_afproj) {
  if (label %in% c('MELD_lambda0','equal')) {
    af <- Reduce('+', lapply(pops, function(p) w[[p]] * bl[[paste0('f_', p)]]))
  } else {
    af <- bl[[paste0('f_', label)]]
  }
  # per-SNP LD score = row-sum of R² (dense within-block)
  ld <- rowSums(R_rep^2)
  data.table(
    chr = bl$chr,
    pos = bl$BP,
    a0  = bl$A2_canon,       # REF
    a1  = bl$A1_canon,       # ALT / effect-coded
    rsid = bl$SNP,
    af   = as.numeric(af),
    ld   = as.numeric(ld),
    group_id = bl$block_id
  )
}

interp_cM <- function(bp, map_bp, map_cM) {
  approx(map_bp, map_cM, xout = bp, rule = 2L)$y
}

# ---------------------------------------------------------------------------
# Per-chromosome loop
map_all <- list()
for (chr in 1:22) {
  chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', chr))
  block_files <- sort(list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE))
  if (length(block_files) == 0L) {
    cat(sprintf('chr%d: no block files, skipping\n', chr)); next
  }

  # Collect per-block R matrices and map rows
  block_R  <- vector('list', length(block_files))
  block_M  <- vector('list', length(block_files))
  for (bi in seq_along(block_files)) {
    bl <- readRDS(block_files[[bi]])
    if (label_needs_all_pops <- (LABEL %in% c('MELD_lambda0','equal'))) {
      if (!all(paste0('R_', pops) %in% names(bl))) {
        stop(sprintf('block %s missing required pops for %s label', block_files[[bi]], LABEL))
      }
    } else {
      if (!(paste0('R_', LABEL) %in% names(bl))) {
        stop(sprintf('block %s missing R_%s for single-pop label', block_files[[bi]], LABEL))
      }
    }
    R <- reconstruct_block(bl, LABEL)
    R_rep <- psd_repair(R)
    block_R[[bi]] <- R_rep
    block_M[[bi]] <- build_map_row(bl, R_rep, LABEL)
  }

  # Stitch to chr-wide block-diagonal sparse matrix
  # Use bdiag from Matrix — it returns a Matrix::CsparseMatrix; forceSymmetric to
  # make it dsCMatrix (symmetric compressed sparse column) which is what
  # bigsnpr::as_SFBM expects.
  R_chr_bd <- bdiag(lapply(block_R, function(m) as(m, 'sparseMatrix')))
  R_chr    <- forceSymmetric(R_chr_bd, uplo = 'U')
  # Cast to dsCMatrix explicitly
  R_chr <- as(R_chr, 'symmetricMatrix')
  R_chr <- as(R_chr, 'CsparseMatrix')
  R_chr <- as(R_chr, 'symmetricMatrix')     # ends up dsCMatrix

  # Save per-chr .rds
  chr_out <- file.path(OUT_DIR, sprintf('LD_with_blocks_chr%d.rds', chr))
  saveRDS(R_chr, chr_out, compress = 'gzip')
  cat(sprintf('chr%d: %d SNPs, %d nonzeros, wrote %s\n',
              chr, nrow(R_chr), length(R_chr@x), chr_out))

  # Combine map rows
  m_chr <- rbindlist(block_M)
  # interpolate cM from ldak_map (LDpred2 doesn't use cM directly, but map format
  # has 'pos' col which is bp — we don't need cM here; LDpred2 only reads chr,pos,a0,a1,af,ld)
  map_all[[chr]] <- m_chr
}

map_dt <- rbindlist(map_all)
map_dt[, af_UKBB := af]    # for backward compat with LDpred2's `names(map)[names(map) == 'af_UKBB']<-'af'`
setcolorder(map_dt, c('chr','pos','a0','a1','rsid','af','af_UKBB','ld','group_id'))
saveRDS(as.data.frame(map_dt), file.path(OUT_DIR, 'map.rds'), compress = 'gzip')
cat(sprintf('map.rds: %d SNPs\n', nrow(map_dt)))
cat('DONE\n')
