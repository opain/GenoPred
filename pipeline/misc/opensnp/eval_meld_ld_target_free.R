#!/usr/bin/env Rscript
# Round 6 Section 4, target-free metrics M1 and M2, evaluated per LDetect block
# on all candidate LD panels vs the analytic meta target R*.
#
# Candidates (per Section 3):
#   EUR, EAS, AFR, CSA, AMR      single-pop empirical (from block RDS)
#   MELD-λ0                     analytic R* — the target; scored for sanity
#   MELD-λ1                     analytic R_pool (mega target, λ=1)
#   MELD-pooled                 empirical LD on pooled EUR+AFR at w=N_p/ΣN_p
#
# Mixture points: w_eur ∈ {0.75, 0.50, 0.25}. Endpoints (w=1, w=0) are omitted
# because R* collapses to R_EUR / R_AFR and every metric is trivial.
#
# Outputs:
#   misc/opensnp/meld_r6_m1.csv          long: w_eur, chr, block_id, m, candidate, err
#   misc/opensnp/meld_r6_m1_summary.csv  per-candidate mean/median/IQR of err
#   misc/opensnp/meld_r6_m2.csv          long: w_eur, chr, block_id, m, candidate, delta, draw, r2
#   misc/opensnp/meld_r6_m2_summary.csv  per-candidate size-weighted mean r² @ each δ

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
refdir       <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
plink2       <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
meld_ld_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld'
sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
out_dir      <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
window_cm    <- 3
w_grid       <- c(0.75, 0.50, 0.25)
N_TOTAL      <- 500000L
n_draws      <- 1L                          # per-block draws; low because n_blocks (~1700) does the averaging
mask_frac    <- 0.10
deltas       <- c(0.01)                     # primary. Add {0.001, 0.1} manually for a sensitivity block later if useful
single_pops  <- c('EUR', 'EAS', 'AFR', 'CSA', 'AMR')
m2_max_m     <- 1200L                       # cap block SNP count for M2 (eigendecomp is n³). M1 uses full block regardless.

args     <- commandArgs(trailingOnly = TRUE)
chr_from <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
chr_to   <- if (length(args) >= 2) as.integer(args[[2]]) else 22L

# ---------------------------------------------------------------------------
# Analytic R at (λ, w) — used for both R* (λ=0) and R_pool (λ=1)
reconstruct_R_lambda <- function(bl, lambda, w_eur, N_total = N_TOTAL) {
  w_afr <- 1 - w_eur
  N_eur <- w_eur * N_total
  N_afr <- w_afr * N_total
  denom <- N_eur + N_afr
  Cov_eur <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_afr <- bl$R_AFR * outer(sqrt(bl$v_AFR), sqrt(bl$v_AFR))
  Sigma_star <- (N_eur * Cov_eur + N_afr * Cov_afr) / denom
  V_star     <- (N_eur * bl$v_EUR + N_afr * bl$v_AFR) / denom
  d <- bl$f_EUR - bl$f_AFR
  B_scalar <- 4 * w_eur * w_afr
  Sigma <- Sigma_star + lambda * B_scalar * outer(d, d)
  V     <- V_star     + lambda * B_scalar * d^2
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

# ---------------------------------------------------------------------------
# Pool composition for MELD-pooled at a given w_eur — maximise pool size while
# preserving the intended n_eur:n_afr ratio, capped by 1KG+HGDP availability.
# For w<0.491 the pool is AFR-capped (take all 688 AFR, subset EUR); for w>0.491
# EUR-capped (take all 665 EUR, subset AFR); at w=0.491 both caps meet.
build_pools <- function(w_grid, seed = 2026L) {
  eur <- readLines(file.path(refdir, 'keep_files', 'EUR.keep'))
  afr <- readLines(file.path(refdir, 'keep_files', 'AFR.keep'))
  set.seed(seed)
  pools <- setNames(vector('list', length(w_grid)),
                    sprintf('w%02d', round(100 * w_grid)))
  for (i in seq_along(w_grid)) {
    w <- w_grid[i]
    if (w <= 0 || w >= 1) {
      # single-pop pool is degenerate; skip pooled candidate at endpoints
      pools[[i]] <- NULL
      next
    }
    ratio_eur_over_afr <- w / (1 - w)
    max_eur_at_ratio   <- length(afr) * ratio_eur_over_afr
    if (max_eur_at_ratio <= length(eur)) {
      # AFR-capped
      n_afr <- length(afr)
      n_eur <- as.integer(round(max_eur_at_ratio))
    } else {
      # EUR-capped
      n_eur <- length(eur)
      n_afr <- as.integer(round(length(eur) / ratio_eur_over_afr))
    }
    eur_pick <- if (n_eur == length(eur)) eur else sample(eur, n_eur)
    afr_pick <- if (n_afr == length(afr)) afr else sample(afr, n_afr)
    pools[[i]] <- list(ids = c(eur_pick, afr_pick),
                       n_eur = length(eur_pick), n_afr = length(afr_pick),
                       w_eur = w,
                       achieved = length(eur_pick) / (length(eur_pick) + length(afr_pick)))
  }
  pools
}

# ---------------------------------------------------------------------------
# Frobenius norm helpers (dense)
frob <- function(M) sqrt(sum(M * M))
frob_err <- function(R_cand, R_target) frob(R_cand - R_target) / frob(R_target)

# M2 core: for a candidate R and vector z, do n_draws mask/impute rounds,
# return a data.table of r² per draw per δ. Reuses the mask across candidates
# so that mask-draw variance doesn't confound candidate comparison.
#
# Efficient version: eigendecompose R_MM once per (candidate, draw). Different
# δ values only offset the eigenvalues, so all δ share the same U/Λ.
run_M2 <- function(R_cand, z, mask_list, deltas, cand_label) {
  out <- vector('list', length(mask_list))
  for (di in seq_along(mask_list)) {
    S <- mask_list[[di]]
    M <- setdiff(seq_along(z), S)
    R_MM <- R_cand[M, M, drop = FALSE]
    R_SM <- R_cand[S, M, drop = FALSE]
    z_M  <- z[M]
    z_S  <- z[S]
    # Eigendecomp; symmetric so this is efficient.
    ev  <- suppressWarnings(eigen(R_MM, symmetric = TRUE))
    U   <- ev$vectors
    lam <- ev$values
    tr_over_n <- sum(lam) / length(lam)   # equal to mean(diag(R_MM))
    Utz <- crossprod(U, z_M)              # UᵀzM, dim = length(M)
    r2s <- numeric(length(deltas))
    for (ki in seq_along(deltas)) {
      ridge <- deltas[[ki]] * tr_over_n
      # z_hat = R_SM @ U @ diag(1/(lam + ridge)) @ UᵀzM
      z_hat <- drop(R_SM %*% (U %*% (Utz / (lam + ridge))))
      r <- suppressWarnings(cor(z_hat, z_S))
      r2s[ki] <- if (is.finite(r)) r^2 else NA_real_
    }
    out[[di]] <- data.table(candidate = cand_label,
                            draw = di, delta = deltas, r2 = r2s)
  }
  rbindlist(out)
}

# ---------------------------------------------------------------------------
# Load per-mixture sumstats z scores once
z_by_w <- list()
for (w in w_grid) {
  tag <- sprintf('eur%02d', round(100 * w))
  ss  <- fread(file.path(sumstats_dir, sprintf('yengo_2022_height_mix_%s.txt', tag)))
  ss[, z := beta / standard_error]
  z_by_w[[sprintf('w%02d', round(100 * w))]] <- ss[, .(SNP = variant_id, z)]
  cat(sprintf('loaded sumstats %s: %d SNPs\n', tag, nrow(ss)))
}

# ---------------------------------------------------------------------------
# Precompute pool sample lists
pools <- build_pools(w_grid)
for (nm in names(pools)) {
  cat(sprintf('pool %s: %d EUR + %d AFR = %d (achieved w_eur = %.4f)\n',
              nm, pools[[nm]]$n_eur, pools[[nm]]$n_afr,
              length(pools[[nm]]$ids), pools[[nm]]$achieved))
}

# ---------------------------------------------------------------------------
# Per-chromosome loop
m1_rows <- list(); m2_rows <- list()
for (chr in chr_from:chr_to) {
  chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', chr))
  block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
  if (length(block_files) == 0L) { cat(sprintf('chr%d: no blocks — skip\n', chr)); next }
  cat(sprintf('\n=== chr%d: %d blocks ===\n', chr, length(block_files)))

  # (a) Attach the chr-wide bed once (for MELD-pooled empirical R)
  tmp <- tempfile(pattern = sprintf('mldeval_chr%d_', chr)); dir.create(tmp)
  bed_prefix <- file.path(tmp, sprintf('chr%d', chr))
  cmd <- c(plink2, '--pfile', file.path(refdir, sprintf('ref.chr%d', chr)),
           '--make-bed', '--out', bed_prefix)
  log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(log, 'status'))) { cat(log, sep='\n'); stop('plink2 failed') }

  bk  <- snp_readBed2(paste0(bed_prefix, '.bed'),
                      backingfile = tempfile(pattern = sprintf('bkev_chr%d_', chr)))
  obj <- snp_attach(bk)
  G   <- obj$genotypes
  fam <- as.data.table(obj$fam)
  bim <- as.data.table(obj$map)
  fam[, row := .I]

  # Per-pool ind.row indices
  ind_row_by_pool <- lapply(pools, function(p) {
    idx <- match(p$ids, fam$sample.ID)
    stopifnot(!any(is.na(idx)))
    idx
  })

  # (b) Per-block loop
  for (fi in seq_along(block_files)) {
    bl <- readRDS(block_files[[fi]])
    m  <- length(bl$SNP)
    if (m < 20L) next                     # skip tiny blocks; masked r² unstable

    # Restrict to SNPs in each per-w sumstats (should overlap heavily — the
    # sumstats were built on the same reference SNP set)
    ord_bim <- match(bl$SNP, bim$marker.ID)
    ok_bim  <- !is.na(ord_bim)
    if (sum(ok_bim) < 20L) next

    # Ridge-align per-w z vector: keep SNPs present in the block and in the
    # per-w sumstats. For simplicity require *all* three w to have the SNPs
    # (they should, since Section 1 wrote the same variant set for w=0.75/0.5/0.25).
    keep_local <- ok_bim
    for (nm in names(z_by_w)) {
      keep_local <- keep_local & bl$SNP %in% z_by_w[[nm]]$SNP
    }
    if (sum(keep_local) < 20L) next
    idx_local <- which(keep_local)
    snp_local <- bl$SNP[idx_local]
    m_use     <- length(idx_local)

    # Subset each of the R matrices to keep_local
    slice_R <- function(R) R[idx_local, idx_local, drop = FALSE]

    # ---- Compute empirical R_pool for each of the 3 pools (using bl$SNP order,
    # then reduce to idx_local afterwards). snp_cor operates in the chr-wide
    # bim indexing (ord_bim); build sign_canon on the same indexing.
    ord_full <- ord_bim
    ord_use  <- ord_full[idx_local]
    a1_bim <- bim$allele1[ord_use]
    a2_bim <- bim$allele2[ord_use]
    alt_e  <- bl$A1_canon[idx_local]
    ref_e  <- bl$A2_canon[idx_local]
    same_al <- (a1_bim == ref_e & a2_bim == alt_e)
    flip_al <- (a1_bim == alt_e & a2_bim == ref_e)
    if (!all(same_al | flip_al)) stop('eval: allele mismatch in slice')
    sign_canon <- ifelse(flip_al, +1, -1)
    pos_M <- bl$cM[idx_local] / 100

    # MELD-pooled empirical R: compute only at the "canonical" w=0.5 pool to bound
    # compute. Full-w sweep for pool-empirical would triple the per-block snp_cor
    # cost and is deferred; the w=0.5 comparison already answers "is empirical R
    # comparable to analytic R_pool at the equal-weight point".
    pool_key_5050 <- 'w50'
    R_pool_emp_50 <- {
      R_sp <- snp_cor(G, ind.row = ind_row_by_pool[[pool_key_5050]],
                      ind.col = ord_use, size = window_cm / 100,
                      infos.pos = pos_M, ncores = 1)
      R <- as.matrix(R_sp)
      diag(R)[is.na(diag(R))] <- 1; R[is.na(R)] <- 0
      R * outer(sign_canon, sign_canon)
    }

    # ---- Prepare z vectors per w
    z_local_by_w <- lapply(names(z_by_w), function(nm) {
      z_dt <- z_by_w[[nm]]
      z_dt[match(snp_local, SNP), z]
    })
    names(z_local_by_w) <- names(z_by_w)

    # ---- Prepare mask draws (shared across candidates and w for identifiability)
    set.seed(2026L + fi)
    mask_list <- lapply(seq_len(n_draws), function(d)
      sort(sample.int(m_use, size = max(2L, round(mask_frac * m_use)))))

    # ---- Loop w_eur values
    for (wi in seq_along(w_grid)) {
      w    <- w_grid[wi]
      wtag <- sprintf('w%02d', round(100 * w))

      R_star <- reconstruct_R_lambda(bl, lambda = 0, w_eur = w)[idx_local, idx_local]
      R_pool <- reconstruct_R_lambda(bl, lambda = 1, w_eur = w)[idx_local, idx_local]

      # Core candidates (always present)
      cand_R <- list(
        MELD_lambda0 = R_star,
        MELD_lambda1 = R_pool,
        EUR          = slice_R(bl$R_EUR),
        AFR          = slice_R(bl$R_AFR)
      )
      # MELD_pooled only at w=0.5 (see comment near R_pool_emp_50)
      if (isTRUE(all.equal(w, 0.5))) cand_R[['MELD_pooled']] <- R_pool_emp_50
      # Optional extra-pop candidates: include per-block only if the block RDS
      # has that pop's R matrix. Blocks without silently drop the candidate.
      for (extra_p in c('EAS', 'CSA', 'AMR')) {
        r_name <- paste0('R_', extra_p)
        if (r_name %in% names(bl)) cand_R[[extra_p]] <- slice_R(bl[[r_name]])
      }

      # M1: Frobenius error to R*
      denom_frob <- frob(R_star)
      for (cn in names(cand_R)) {
        num <- frob(cand_R[[cn]] - R_star)
        err <- if (denom_frob > 0) num / denom_frob else NA_real_
        m1_rows[[length(m1_rows) + 1L]] <- data.table(
          w_eur = w, chr = chr, block_id = bl$block_id,
          m = m_use, candidate = cn, err = err
        )
      }

      # M2: masked-z re-imputation
      # Skip MELD-λ0: R* is the target, so masked-z imputation using R* is
      # a self-consistency check. Included in M1 for completeness (err=0 by
      # construction), but M2 near-1 is trivial and burns compute.
      z_here <- z_local_by_w[[wtag]]
      if (any(!is.finite(z_here))) next   # NA z means missing SNP from sumstats — skip

      # Cap block SNP count for M2 to bound eigendecomp cost (n^3). Deterministic
      # subsample per block. Uses the same subset across all candidates and w so
      # comparisons stay apples-to-apples within the block.
      if (m_use > m2_max_m) {
        set.seed(2027L + fi)
        m2_idx <- sort(sample.int(m_use, m2_max_m))
        cand_R_m2 <- lapply(cand_R, function(R) R[m2_idx, m2_idx, drop = FALSE])
        z_m2      <- z_here[m2_idx]
        m_use_m2  <- m2_max_m
      } else {
        cand_R_m2 <- cand_R
        z_m2      <- z_here
        m_use_m2  <- m_use
      }
      # Regenerate masks scoped to (possibly-subsetted) m_use_m2
      set.seed(2028L + fi)
      mask_list_m2 <- lapply(seq_len(n_draws), function(d)
        sort(sample.int(m_use_m2, size = max(2L, round(mask_frac * m_use_m2)))))

      for (cn in setdiff(names(cand_R_m2), 'MELD_lambda0')) {
        rows <- run_M2(cand_R_m2[[cn]], z_m2, mask_list_m2, deltas, cn)
        rows[, `:=`(w_eur = w, chr = chr, block_id = bl$block_id, m = m_use_m2)]
        m2_rows[[length(m2_rows) + 1L]] <- rows
      }
    }

    if (fi %% 20L == 0L) cat(sprintf('  chr%d: %d/%d blocks done\n', chr, fi, length(block_files)))
  }

  # cleanup chr backing files
  bk_files <- c(bk, sub('\\.bk$', '.rds', bk))
  suppressWarnings(file.remove(bk_files[file.exists(bk_files)]))
  unlink(tmp, recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Write outputs
m1 <- rbindlist(m1_rows)
m2 <- rbindlist(m2_rows)
fwrite(m1, file.path(out_dir, 'meld_r6_m1.csv'))
fwrite(m2, file.path(out_dir, 'meld_r6_m2.csv'))
cat(sprintf('wrote %s (%d rows)\n', file.path(out_dir, 'meld_r6_m1.csv'), nrow(m1)))
cat(sprintf('wrote %s (%d rows)\n', file.path(out_dir, 'meld_r6_m2.csv'), nrow(m2)))

# M1 summary: per (w_eur, candidate)
m1_summary <- m1[, .(n_blocks = .N,
                     err_mean   = mean(err,   na.rm = TRUE),
                     err_median = median(err, na.rm = TRUE),
                     err_q25    = quantile(err, .25, na.rm = TRUE),
                     err_q75    = quantile(err, .75, na.rm = TRUE)),
                 by = .(w_eur, candidate)]
fwrite(m1_summary, file.path(out_dir, 'meld_r6_m1_summary.csv'))

# M2 summary: mean r² per (w_eur, candidate, delta), size-weighted by m
m2_summary <- m2[, .(n_blocks = length(unique(paste(chr, block_id))),
                     n_draws  = .N,
                     r2_mean_unweighted = mean(r2, na.rm = TRUE),
                     r2_mean_wm        = sum(r2 * m, na.rm = TRUE) / sum(m * (!is.na(r2)))),
                 by = .(w_eur, candidate, delta)]
fwrite(m2_summary, file.path(out_dir, 'meld_r6_m2_summary.csv'))

cat('\n=== M1 summary (mean err by w & candidate) ===\n')
print(dcast(m1_summary, candidate ~ w_eur, value.var = 'err_mean'))
cat('\n=== M2 summary at delta = 0.01 (size-weighted mean r² by w & candidate) ===\n')
m2_01 <- m2_summary[delta == 0.01]
print(dcast(m2_01, candidate ~ w_eur, value.var = 'r2_mean_wm'))

cat('\nDONE\n')
