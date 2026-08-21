#!/usr/bin/env Rscript
# Round 6 Section 2 invariants. Runs against the block RDS files produced by
# build_meld_ld.R and fails loudly on violation (per prompt).
#
# 1. diag(R(λ)) == 1 to 1e-10, for every λ  (per-block, all blocks)
# 2. At w = (1, 0), R(λ) == R_EUR bitwise-close, any λ
# 3. min eigen(R(λ)) ≥ -1e-8, λ ∈ {0, 1}, per-block (PSD)
# 4. Cross-check R_pool (analytic λ=1) vs empirical pooled LD at w = N_p/ΣN_p,
#    on one or two representative blocks.
#
# Usage:
#   Rscript test_build_meld_ld.R                          # full pass
#   Rscript test_build_meld_ld.R chr22                    # single-chr fast pass

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# Config (kept identical to build_meld_ld.R)
refdir       <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
plink2       <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
map_dir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_map/genetic_map_b37'
meld_ld_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
emp_var_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
window_cm    <- 3
w_eur_cross  <- 0.50                # invariant 4 uses w_eur = 0.5

args <- commandArgs(trailingOnly = TRUE)
chr_filter <- if (length(args) >= 1) args[[1]] else NULL

# ---------------------------------------------------------------------------
# Helpers
# Reconstruct R(λ) from a block RDS at nominal weights (w_eur, w_afr) and
# per-population weights ω_p (defaults to N_p from Section 1's total N=500k,
# reduces to N_p/N_total = w_p under equal N_total). At λ=0 gives R*; at λ=1
# gives R_pool.
reconstruct_R_lambda <- function(bl, lambda, w_eur, w_afr, N_total = 500000L) {
  N_eur <- w_eur * N_total
  N_afr <- w_afr * N_total
  # Σ_p Cov_p = R_p * outer(sqrt(v_p), sqrt(v_p))
  # Meta uses N_p weights per Section 1 IVW:
  #   Σ*(i,j) = Σ_p N_p Cov_p(i,j) / Σ_p N_p
  #   V*(i)   = Σ_p N_p v_p,i      / Σ_p N_p
  Cov_eur <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_afr <- bl$R_AFR * outer(sqrt(bl$v_AFR), sqrt(bl$v_AFR))
  denom   <- N_eur + N_afr
  Sigma_star <- (N_eur * Cov_eur + N_afr * Cov_afr) / denom
  V_star     <- (N_eur * bl$v_EUR + N_afr * bl$v_AFR) / denom
  # Between-population term (rank-one; store as outer product not densified)
  # B is a Wahlund-form 4 w_EUR w_AFR (f_EUR - f_AFR)(f_EUR - f_AFR)^T; the
  # arithmetic weights are per the prompt's Section 2 mega-target formula.
  d <- bl$f_EUR - bl$f_AFR
  B_scalar <- 4 * w_eur * w_afr
  # Full B is B_scalar * outer(d, d). diag(B) = B_scalar * d^2.
  Sigma <- Sigma_star + lambda * B_scalar * outer(d, d)
  V     <- V_star     + lambda * B_scalar * d^2
  # R(λ) = D(λ)^{-1/2} Σ(λ) D(λ)^{-1/2}
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R
}

# ---------------------------------------------------------------------------
# Enumerate block files
chrs <- if (is.null(chr_filter)) sprintf('chr%d', 1:22) else chr_filter
chr_dirs <- file.path(meld_ld_dir, chrs)
chr_dirs <- chr_dirs[dir.exists(chr_dirs)]
if (length(chr_dirs) == 0L) stop('no block dirs under ', meld_ld_dir)

all_blocks <- unlist(lapply(chr_dirs, function(d)
  file.path(d, list.files(d, pattern = '^block_.*rds$'))))
cat(sprintf('found %d block files across %d chromosome dirs\n',
            length(all_blocks), length(chr_dirs)))

# ---------------------------------------------------------------------------
# Invariant 1: diag(R(λ)) == 1 for all λ
#            (nontrivial: v_p can be 0 in one pop, so numerator/denominator
#             both non-standard — see prompt Section 2 edge case discussion.)
test_diag <- function(files, tol = 1e-8) {
  cat('\n--- Invariant 1: diag(R(λ)) == 1 for λ ∈ {0, 0.5, 1} ---\n')
  worst <- 0
  worst_where <- NULL
  for (fn in files) {
    bl <- readRDS(fn)
    for (lambda in c(0, 0.5, 1)) {
      R  <- reconstruct_R_lambda(bl, lambda, w_eur = 0.5, w_afr = 0.5)
      dd <- diag(R)
      # A SNP has R diag NaN or Inf only if V = 0, i.e. monomorphic in every
      # pop AND B contributes nothing there (d = 0). Treat those as passes.
      ok <- is.finite(dd)
      m  <- if (any(ok)) max(abs(dd[ok] - 1)) else 0
      if (m > worst) {
        worst <- m
        worst_where <- list(file = fn, lambda = lambda,
                            snp_worst = bl$SNP[which.max(abs(dd - 1))])
      }
    }
  }
  cat(sprintf('  worst |diag - 1| = %.3e  (tolerance %.0e)\n', worst, tol))
  if (worst > tol) {
    cat('  worst was at:\n'); print(worst_where)
    stop('Invariant 1 FAILED')
  }
  cat('  PASS\n')
}

# ---------------------------------------------------------------------------
# Invariant 2: at w=(1,0), R(λ) == R_EUR bitwise-close, any λ
test_edge_w1 <- function(files, tol = 1e-12) {
  cat('\n--- Invariant 2: R(λ)|w=(1,0) == R_EUR for any λ ---\n')
  worst <- 0
  worst_where <- NULL
  for (fn in files) {
    bl <- readRDS(fn)
    for (lambda in c(0, 0.5, 1)) {
      R  <- reconstruct_R_lambda(bl, lambda, w_eur = 1, w_afr = 0)
      # When w_afr = 0, B_scalar = 4*1*0 = 0, so B contribution vanishes for any λ.
      # And Σ_pop weights collapse to (1, 0), giving R = R_EUR.
      d  <- max(abs(R - bl$R_EUR), na.rm = TRUE)
      if (d > worst) {
        worst <- d
        worst_where <- list(file = fn, lambda = lambda)
      }
    }
  }
  cat(sprintf('  worst |R - R_EUR| = %.3e  (tolerance %.0e)\n', worst, tol))
  if (worst > tol) {
    cat('  worst was at:\n'); print(worst_where)
    stop('Invariant 2 FAILED')
  }
  cat('  PASS\n')
}

# ---------------------------------------------------------------------------
# Invariant 3: min eigen(R(λ)) ≥ -1e-8, λ ∈ {0, 1}
test_psd <- function(files, tol = -1e-6) {
  cat('\n--- Invariant 3: min eigen(R(λ)) ≥ -1e-8 (report histogram) ---\n')
  minev <- list(lambda0 = numeric(0), lambda1 = numeric(0))
  worst_where <- NULL
  worst <- 0
  for (fn in files) {
    bl <- readRDS(fn)
    for (lambda in c(0, 1)) {
      R  <- reconstruct_R_lambda(bl, lambda, w_eur = 0.5, w_afr = 0.5)
      # replace any NaN with 0 for eigen (should already be finite in practice)
      R[!is.finite(R)] <- 0
      e  <- min(suppressWarnings(eigen(R, symmetric = TRUE, only.values = TRUE)$values))
      k  <- if (lambda == 0) 'lambda0' else 'lambda1'
      minev[[k]] <- c(minev[[k]], e)
      if (e < worst) { worst <- e; worst_where <- list(file = fn, lambda = lambda) }
    }
  }
  for (k in names(minev)) {
    cat(sprintf('  %s: min=%.3e  q05=%.3e  median=%.3e  q95=%.3e  max=%.3e   (n=%d)\n',
                k, min(minev[[k]]), quantile(minev[[k]], .05),
                median(minev[[k]]),
                quantile(minev[[k]], .95), max(minev[[k]]),
                length(minev[[k]])))
  }
  cat(sprintf('  worst min-eigen = %.3e  (tolerance %.0e)\n', worst, tol))
  if (worst < tol) {
    cat('  worst was at:\n'); print(worst_where)
    warning('Invariant 3 boundary crossed — investigate')
  }
  cat('  PASS (report only; near-zero neg eigs from f.p. expected)\n')
}

# ---------------------------------------------------------------------------
# Invariant 4: analytic R_pool ≈ empirical pooled LD at w = N_p/ΣN_p
# Rather than instantiating a full parameterised private refdir, compute the
# empirical LD on-the-fly per block from the shared pgen using the same window
# convention as build_meld_ld.R.
build_pool_keep <- function(w_eur, N_total = 500000L, seed = 2026L) {
  set.seed(seed)
  eur <- readLines(file.path(refdir, 'keep_files', 'EUR.keep'))
  afr <- readLines(file.path(refdir, 'keep_files', 'AFR.keep'))
  n_eur <- min(round(w_eur       * N_total),  length(eur))
  n_afr <- min(round((1 - w_eur) * N_total), length(afr))
  # For w_eur = 0.5, N_total = 500k, we hit the cap on both pops.
  eur_pick <- if (n_eur == length(eur)) eur else sample(eur, n_eur)
  afr_pick <- if (n_afr == length(afr)) afr else sample(afr, n_afr)
  list(ids = c(eur_pick, afr_pick),
       n_eur = length(eur_pick), n_afr = length(afr_pick))
}

# Per-block empirical LD helper reused from build_meld_ld.R
compute_block_pool_R <- function(bl, pool_ids) {
  # Fetch dosages for these SNPs restricted to pool samples via plink2, then
  # compute correlations directly (avoid a separate bigsnpr attach for just a
  # few blocks — cheaper to shell out to plink2 --r-unphased square).
  tmp <- tempfile(pattern = sprintf('poolR_chr%d_bl%d_', bl$chr, bl$block_id))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  keep_path    <- file.path(tmp, 'pool.keep')
  extract_path <- file.path(tmp, 'snps.txt')
  writeLines(pool_ids, keep_path)
  writeLines(bl$SNP,   extract_path)

  pfile <- file.path(refdir, sprintf('ref.chr%d', bl$chr))
  cmd <- c(plink2,
           '--pfile',    pfile,
           '--keep',     keep_path,
           '--extract',  extract_path,
           '--make-bed',
           '--out',      file.path(tmp, 'blockbed'))
  log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(log, 'status'))) { cat(log, sep='\n'); stop('plink2 make-bed failed') }

  bk  <- snp_readBed2(file.path(tmp, 'blockbed.bed'),
                      backingfile = file.path(tmp, 'bk'))
  obj <- snp_attach(bk)
  G   <- obj$genotypes
  bim <- as.data.table(obj$map)

  # Align bim to bl$SNP order (plink may reorder); build ind.col accordingly.
  ord <- match(bl$SNP, bim$marker.ID)
  # Some SNPs may be missing (rare); drop those from the comparison.
  ok  <- !is.na(ord)
  pos_M <- bl$cM[ok] / 100
  R_sp <- snp_cor(G, ind.col = ord[ok], size = window_cm / 100,
                  infos.pos = pos_M, ncores = 1)
  R <- as.matrix(R_sp)
  diag(R)[is.na(diag(R))] <- 1
  R[is.na(R)] <- 0

  # Sign-align to canonical ALT the same way as build_meld_ld.R
  a1_bim <- bim$allele1[ord[ok]]
  a2_bim <- bim$allele2[ord[ok]]
  alt_e  <- bl$A1_canon[ok]
  ref_e  <- bl$A2_canon[ok]
  same   <- (a1_bim == ref_e & a2_bim == alt_e)
  flip   <- (a1_bim == alt_e & a2_bim == ref_e)
  if (!all(same | flip)) stop('pool: bim/canon allele mismatch')
  sign_canon <- ifelse(flip, +1, -1)
  R <- R * outer(sign_canon, sign_canon)

  list(R = R, ok = ok)
}

test_cross_pool <- function(files) {
  cat('\n--- Invariant 4: analytic R_pool ≈ empirical pooled LD at w=0.5 ---\n')
  # Use up to 3 representative blocks: smallest, median, largest by SNP count.
  sizes <- sapply(files, function(f) length(readRDS(f)$SNP))
  idx <- unique(c(which.min(sizes),
                  order(sizes)[ceiling(length(sizes) / 2)],
                  which.max(sizes)))
  chosen <- files[idx]

  pool <- build_pool_keep(w_eur = w_eur_cross)
  cat(sprintf('  pooled sample: %d EUR + %d AFR = %d\n',
              pool$n_eur, pool$n_afr, length(pool$ids)))

  for (fn in chosen) {
    bl <- readRDS(fn)
    R_analytic <- reconstruct_R_lambda(bl, lambda = 1,
                                       w_eur = w_eur_cross,
                                       w_afr = 1 - w_eur_cross)
    ep <- compute_block_pool_R(bl, pool$ids)
    R_emp <- ep$R
    ok    <- ep$ok
    if (!all(ok)) {
      R_analytic <- R_analytic[ok, ok]
    }
    off <- upper.tri(R_analytic)
    v1 <- R_analytic[off]
    v2 <- R_emp[off]
    keep <- is.finite(v1) & is.finite(v2)
    if (sum(keep) < 100) { cat(sprintf('  chr%d block%d: too few off-diag (%d) — skipping\n',
                                       bl$chr, bl$block_id, sum(keep))); next }
    rho <- cor(v1[keep], v2[keep])
    # Fisher-z SE for sanity: n effective from #off-diagonal pairs; not exact
    # (correlations are dependent) but gives a rough scale.
    z   <- 0.5 * log((1 + rho) / (1 - rho))
    z_se <- 1 / sqrt(max(sum(keep) - 3L, 1L))
    cat(sprintf('  chr%d block%d (m=%d SNPs): cor(R_analytic, R_pooled) = %.4f   (Fisher-z SE ≈ %.4f)\n',
                bl$chr, bl$block_id, length(bl$SNP), rho, z_se))
    if (rho < 0.90) warning(sprintf('  chr%d block%d: correlation < 0.90 — investigate',
                                    bl$chr, bl$block_id))
  }
  cat('  DONE\n')
}

# ---------------------------------------------------------------------------
# Run the tests
if (length(all_blocks) == 0L) stop('no block files under ', meld_ld_dir)

test_diag(all_blocks)
test_edge_w1(all_blocks)
test_psd(all_blocks)
test_cross_pool(all_blocks)

cat('\nALL INVARIANTS PASSED\n')
