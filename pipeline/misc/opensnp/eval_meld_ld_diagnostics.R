#!/usr/bin/env Rscript
# Round 6b diagnostics — chr22 only.
#
# Investigates whether the checkpoint's M2 finding ("MELD-λ1 beats every
# single-population candidate") is a real ancestry effect or an artefact of
#   (i)  omitting MELD-λ0 from the M2 comparison
#   (ii) conditioning: windowed snp_cor is indefinite; B is a positive rank-one
#        term, so adding it repairs conditioning
#   (iii) circularity: SE_meta and B share the reference AF ingredients.
#
# Diagnostics implemented:
#   D1  Enter MELD-λ0 in M2 (no exclusions)
#   D2  Uniform PSD repair (eigen-clamp + diag-rescale) for every candidate;
#       report pre/post min eigenvalue, pre/post M2, and rank correlation of
#       pre-repair min-eigen vs pre-repair M2 across candidates
#   D3  Placebo — permute ΔAF vector within block, form rank-one term the same
#       way B is formed; several draws
#   D4  λ sweep {0, 0.25, 0.5, 0.75, 1} post-repair
#   D5  Stratify (M2_λ=1 − M2_λ=0) by |ΔAF| tertile within each block
#   D6  handled post-hoc in a separate bootstrap script.
#
# Writes ALL outputs to new files (do not overwrite meld_r6_m1*.csv / meld_r6_m2*.csv):
#   meld_r6b_d1_m2_all_candidates.csv
#   meld_r6b_d2_repair.csv
#   meld_r6b_d3_placebo.csv
#   meld_r6b_d4_lambda_sweep.csv
#   meld_r6b_d5_af_stratified.csv
#
# Usage: Rscript eval_meld_ld_diagnostics.R [chr]         # default: 22

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
lambda_grid  <- c(0, 0.25, 0.5, 0.75, 1)   # D4
deltas       <- c(0.001, 0.01, 0.1)        # for D2 sensitivity + D1 primary at 0.01
n_perm       <- 8L                         # D3 placebo draws
n_draws_mask <- 1L                         # per-block mask draws (blocks average noise; 1 was Round-6 baseline)
mask_frac    <- 0.10
N_TOTAL      <- 500000L                    # matches Section 1
n_af_bins    <- 3L                         # D5 tertile stratification

args <- commandArgs(trailingOnly = TRUE)
CHR  <- if (length(args) >= 1) as.integer(args[[1]]) else 22L

# ---------------------------------------------------------------------------
# Reconstruction identical to Section 2 build_meld_ld / eval_meld_ld_target_free
reconstruct_R_lambda <- function(bl, lambda, w_eur, N_total = N_TOTAL,
                                 d_override = NULL) {
  w_afr <- 1 - w_eur
  N_eur <- w_eur * N_total
  N_afr <- w_afr * N_total
  denom <- N_eur + N_afr
  Cov_eur <- bl$R_EUR * outer(sqrt(bl$v_EUR), sqrt(bl$v_EUR))
  Cov_afr <- bl$R_AFR * outer(sqrt(bl$v_AFR), sqrt(bl$v_AFR))
  Sigma_star <- (N_eur * Cov_eur + N_afr * Cov_afr) / denom
  V_star     <- (N_eur * bl$v_EUR + N_afr * bl$v_AFR) / denom
  d <- if (is.null(d_override)) bl$f_EUR - bl$f_AFR else d_override
  B_scalar <- 4 * w_eur * w_afr
  Sigma <- Sigma_star + lambda * B_scalar * outer(d, d)
  V     <- V_star     + lambda * B_scalar * d^2
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

# PSD repair: eigen-clamp to zero, reconstruct, rescale diagonal to 1
psd_repair <- function(R) {
  ev <- eigen(R, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  R2 <- ev$vectors %*% (lam * t(ev$vectors))
  # cov2cor-style rescale
  d <- sqrt(pmax(diag(R2), 0))
  d[d == 0] <- 1
  R2 <- R2 / outer(d, d)
  diag(R2) <- 1
  R2
}

# Min eigenvalue helper (only values, not vectors — faster)
min_eig <- function(R) {
  suppressWarnings(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values))
}

# ---------------------------------------------------------------------------
# M2 core reused: eigendecomp of R_MM (used across the deltas we sweep)
run_M2_one <- function(R_cand, z, mask_S, deltas) {
  m_all <- seq_along(z)
  M <- setdiff(m_all, mask_S)
  R_MM <- R_cand[M, M, drop = FALSE]
  R_SM <- R_cand[mask_S, M, drop = FALSE]
  z_M  <- z[M]; z_S  <- z[mask_S]
  ev  <- suppressWarnings(eigen(R_MM, symmetric = TRUE))
  U   <- ev$vectors; lam <- ev$values
  tr_over_n <- sum(lam) / length(lam)
  Utz <- crossprod(U, z_M)
  sapply(deltas, function(dlt) {
    ridge <- dlt * abs(tr_over_n)      # abs to be robust to sign-flipped tr
    z_hat <- drop(R_SM %*% (U %*% (Utz / (lam + ridge))))
    r <- suppressWarnings(cor(z_hat, z_S))
    if (is.finite(r)) r^2 else NA_real_
  })
}

# ---------------------------------------------------------------------------
# Load per-mixture sumstats z-scores
z_by_w <- list()
for (w in w_grid) {
  tag <- sprintf('eur%02d', round(100 * w))
  ss  <- fread(file.path(sumstats_dir, sprintf('yengo_2022_height_mix_%s.txt', tag)))
  ss[, z := beta / standard_error]
  z_by_w[[sprintf('w%02d', round(100 * w))]] <- ss[, .(SNP = variant_id, z)]
}
cat(sprintf('loaded sumstats for %d w values\n', length(z_by_w)))

# ---------------------------------------------------------------------------
# Empirical MELD-pooled at w=0.5 (chr-wide bed attach; per-block ind.col slice)
# Also used for the pool ind.row across the block loop.
build_pool_w50 <- function(seed = 2026L) {
  set.seed(seed)
  eur <- readLines(file.path(refdir, 'keep_files', 'EUR.keep'))
  afr <- readLines(file.path(refdir, 'keep_files', 'AFR.keep'))
  n_eur <- min(length(eur), length(afr))
  n_afr <- n_eur
  eur_pick <- if (n_eur == length(eur)) eur else sample(eur, n_eur)
  afr_pick <- if (n_afr == length(afr)) afr else sample(afr, n_afr)
  list(ids = c(eur_pick, afr_pick), n_eur = length(eur_pick), n_afr = length(afr_pick))
}
pool_50 <- build_pool_w50()
cat(sprintf('pool w50: %d EUR + %d AFR\n', pool_50$n_eur, pool_50$n_afr))

# Attach chr bed once
tmp <- tempfile(pattern = sprintf('mldiag_chr%d_', CHR)); dir.create(tmp)
bed_prefix <- file.path(tmp, sprintf('chr%d', CHR))
cmd <- c(plink2, '--pfile', file.path(refdir, sprintf('ref.chr%d', CHR)),
         '--make-bed', '--out', bed_prefix)
log <- system2(cmd[1], cmd[-1], stdout = TRUE, stderr = TRUE)
if (!is.null(attr(log, 'status'))) { cat(log, sep='\n'); stop('plink2 failed') }
bk  <- snp_readBed2(paste0(bed_prefix, '.bed'),
                    backingfile = tempfile(pattern = sprintf('bkd_chr%d_', CHR)))
obj <- snp_attach(bk)
G <- obj$genotypes; fam <- as.data.table(obj$fam); bim <- as.data.table(obj$map)
fam[, row := .I]
ind_row_pool <- match(pool_50$ids, fam$sample.ID)
stopifnot(!any(is.na(ind_row_pool)))

# ---------------------------------------------------------------------------
# Per-block loop
chr_dir <- file.path(meld_ld_dir, sprintf('chr%d', CHR))
block_files <- list.files(chr_dir, pattern = '^block_.*rds$', full.names = TRUE)
cat(sprintf('chr%d: %d blocks\n', CHR, length(block_files)))

d1_rows <- list(); d2_rows <- list(); d3_rows <- list()
d4_rows <- list(); d5_rows <- list()

for (fi in seq_along(block_files)) {
  bl <- readRDS(block_files[[fi]])
  m_full <- length(bl$SNP)

  # Restrict to SNPs also in all per-w sumstats
  keep <- rep(TRUE, m_full)
  for (nm in names(z_by_w)) keep <- keep & bl$SNP %in% z_by_w[[nm]]$SNP
  if (sum(keep) < 30L) next
  idx <- which(keep)

  ord_bim <- match(bl$SNP, bim$marker.ID)
  ok_bim  <- !is.na(ord_bim); idx <- intersect(idx, which(ok_bim))
  if (length(idx) < 30L) next

  # Alignment for pool empirical R (same as build_meld_ld.R canonical ALT)
  a1_bim <- bim$allele1[ord_bim[idx]]; a2_bim <- bim$allele2[ord_bim[idx]]
  alt_e  <- bl$A1_canon[idx];          ref_e  <- bl$A2_canon[idx]
  same_al <- (a1_bim == ref_e & a2_bim == alt_e)
  flip_al <- (a1_bim == alt_e & a2_bim == ref_e)
  if (!all(same_al | flip_al)) stop('allele mismatch in slice')
  sign_canon <- ifelse(flip_al, +1, -1)
  pos_M <- bl$cM[idx] / 100

  # Slice block ingredients to idx
  slice_bl <- list(
    chr = bl$chr, block_id = bl$block_id, SNP = bl$SNP[idx],
    A1_canon = bl$A1_canon[idx], A2_canon = bl$A2_canon[idx],
    v_EUR = bl$v_EUR[idx], v_AFR = bl$v_AFR[idx],
    f_EUR = bl$f_EUR[idx], f_AFR = bl$f_AFR[idx],
    R_EUR = bl$R_EUR[idx, idx, drop = FALSE],
    R_AFR = bl$R_AFR[idx, idx, drop = FALSE]
  )

  # Determine mask (shared across candidates/w for identifiability)
  m_use <- length(idx)
  set.seed(2028L + fi)
  mask_S <- sort(sample.int(m_use, size = max(2L, round(mask_frac * m_use))))

  # z vectors per w
  z_by_w_local <- lapply(names(z_by_w), function(nm) {
    z_by_w[[nm]][match(slice_bl$SNP, SNP), z]
  })
  names(z_by_w_local) <- names(z_by_w)
  if (any(sapply(z_by_w_local, function(z) any(!is.finite(z))))) next

  # Compute R_pool_empirical at w=0.5 for this block
  R_pool_emp <- {
    R_sp <- snp_cor(G, ind.row = ind_row_pool, ind.col = ord_bim[idx],
                    size = window_cm / 100, infos.pos = pos_M, ncores = 1)
    R <- as.matrix(R_sp)
    diag(R)[is.na(diag(R))] <- 1; R[is.na(R)] <- 0
    R * outer(sign_canon, sign_canon)
  }

  # ---- Build the candidate R matrices at each w (for D1, D2) -----------
  for (wi in seq_along(w_grid)) {
    w <- w_grid[wi]
    R_lambda0 <- reconstruct_R_lambda(slice_bl, 0, w)
    R_lambda1 <- reconstruct_R_lambda(slice_bl, 1, w)
    cand_R_raw <- list(
      MELD_lambda0 = R_lambda0,
      MELD_lambda1 = R_lambda1,
      EUR          = slice_bl$R_EUR,
      AFR          = slice_bl$R_AFR
    )
    if (isTRUE(all.equal(w, 0.5))) cand_R_raw$MELD_pooled <- R_pool_emp

    z_here <- z_by_w_local[[sprintf('w%02d', round(100 * w))]]

    for (cn in names(cand_R_raw)) {
      R <- cand_R_raw[[cn]]
      me_pre  <- min_eig(R)
      # D1 M2 pre-repair
      r2_pre  <- run_M2_one(R, z_here, mask_S, deltas)

      # D2 PSD repair
      R_rep <- psd_repair(R)
      me_post <- min_eig(R_rep)
      r2_post <- run_M2_one(R_rep, z_here, mask_S, deltas)

      d1_rows[[length(d1_rows) + 1L]] <- data.table(
        chr = slice_bl$chr, block_id = slice_bl$block_id, m = m_use,
        w_eur = w, candidate = cn,
        delta = deltas, r2_pre = r2_pre, r2_post = r2_post
      )
      d2_rows[[length(d2_rows) + 1L]] <- data.table(
        chr = slice_bl$chr, block_id = slice_bl$block_id, m = m_use,
        w_eur = w, candidate = cn,
        min_eig_pre = me_pre, min_eig_post = me_post,
        r2_pre = r2_pre[which(deltas == 0.01)],
        r2_post = r2_post[which(deltas == 0.01)]
      )
    }
  }

  # ---- D3: placebo — permute ΔAF, form the rank-one term, score M2 -----
  d_true <- slice_bl$f_EUR - slice_bl$f_AFR
  for (perm_i in seq_len(n_perm)) {
    set.seed(3000L + fi * 100L + perm_i)
    d_perm <- sample(d_true)
    for (w in w_grid) {
      # Reconstruct λ=1 with the permuted d (Wahlund weight and formula unchanged)
      R_placebo <- reconstruct_R_lambda(slice_bl, 1, w, d_override = d_perm)
      R_placebo <- psd_repair(R_placebo)   # apply repair uniformly, per D2
      z_here <- z_by_w_local[[sprintf('w%02d', round(100 * w))]]
      r2 <- run_M2_one(R_placebo, z_here, mask_S, deltas)
      d3_rows[[length(d3_rows) + 1L]] <- data.table(
        chr = slice_bl$chr, block_id = slice_bl$block_id, m = m_use,
        w_eur = w, perm = perm_i, delta = deltas, r2 = r2
      )
    }
  }

  # ---- D4: λ sweep, post-repair ---------------------------------------
  for (w in w_grid) {
    for (lam in lambda_grid) {
      R_lam <- reconstruct_R_lambda(slice_bl, lam, w)
      R_lam <- psd_repair(R_lam)
      z_here <- z_by_w_local[[sprintf('w%02d', round(100 * w))]]
      r2 <- run_M2_one(R_lam, z_here, mask_S, deltas)
      d4_rows[[length(d4_rows) + 1L]] <- data.table(
        chr = slice_bl$chr, block_id = slice_bl$block_id, m = m_use,
        w_eur = w, lambda = lam, delta = deltas, r2 = r2
      )
    }
  }

  # ---- D5: stratify (M2_λ=1 − M2_λ=0) by |ΔAF| tertile within block ----
  # Bin SNPs by |ΔAF|, then compute M2 restricted to each bin as *both* the
  # kept and masked set. Approach: within each bin, redraw mask restricted to
  # bin SNPs; use only bin SNPs' rows/cols of R for the M2 solve.
  abs_d <- abs(d_true)
  bin_cuts <- quantile(abs_d, probs = seq(0, 1, length.out = n_af_bins + 1L))
  bin_id <- as.integer(cut(abs_d, breaks = bin_cuts, include.lowest = TRUE))
  bin_id[is.na(bin_id)] <- 1L
  for (w in w_grid) {
    R_l0 <- psd_repair(reconstruct_R_lambda(slice_bl, 0, w))
    R_l1 <- psd_repair(reconstruct_R_lambda(slice_bl, 1, w))
    z_here <- z_by_w_local[[sprintf('w%02d', round(100 * w))]]
    for (b in seq_len(n_af_bins)) {
      idx_b <- which(bin_id == b)
      if (length(idx_b) < 30L) next
      # Restricted mask within the bin
      set.seed(5000L + fi * 100L + b)
      mask_b <- sort(sample.int(length(idx_b),
                                size = max(2L, round(mask_frac * length(idx_b)))))
      r2_l0 <- run_M2_one(R_l0[idx_b, idx_b, drop = FALSE], z_here[idx_b],
                          mask_b, deltas)
      r2_l1 <- run_M2_one(R_l1[idx_b, idx_b, drop = FALSE], z_here[idx_b],
                          mask_b, deltas)
      d5_rows[[length(d5_rows) + 1L]] <- data.table(
        chr = slice_bl$chr, block_id = slice_bl$block_id,
        w_eur = w, af_bin = b, m_bin = length(idx_b),
        af_bin_lo = bin_cuts[b], af_bin_hi = bin_cuts[b + 1L],
        delta = deltas, r2_l0 = r2_l0, r2_l1 = r2_l1, diff = r2_l1 - r2_l0
      )
    }
  }

  if (fi %% 5L == 0L) cat(sprintf('  chr%d: %d/%d blocks done\n', CHR, fi, length(block_files)))
}

# ---------------------------------------------------------------------------
# Write outputs
d1 <- rbindlist(d1_rows); d2 <- rbindlist(d2_rows); d3 <- rbindlist(d3_rows)
d4 <- rbindlist(d4_rows); d5 <- rbindlist(d5_rows)

fwrite(d1, file.path(out_dir, 'meld_r6b_d1_m2_all_candidates.csv'))
fwrite(d2, file.path(out_dir, 'meld_r6b_d2_repair.csv'))
fwrite(d3, file.path(out_dir, 'meld_r6b_d3_placebo.csv'))
fwrite(d4, file.path(out_dir, 'meld_r6b_d4_lambda_sweep.csv'))
fwrite(d5, file.path(out_dir, 'meld_r6b_d5_af_stratified.csv'))

cat(sprintf('\nwrote D1 (%d rows) D2 (%d) D3 (%d) D4 (%d) D5 (%d)\n',
            nrow(d1), nrow(d2), nrow(d3), nrow(d4), nrow(d5)))

# ---------------------------------------------------------------------------
# Quick console summaries (full analysis + bootstrap in a separate script)
cat('\n=== D1: M2 by candidate × w (size-weighted mean r², δ=0.01, PRE-repair) ===\n')
sw <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w * !is.na(x))
d1_01 <- d1[delta == 0.01]
tab_pre  <- d1_01[, .(r2_pre_wm  = sw(r2_pre,  m)), by = .(candidate, w_eur)]
tab_post <- d1_01[, .(r2_post_wm = sw(r2_post, m)), by = .(candidate, w_eur)]
print(dcast(tab_pre,  candidate ~ w_eur, value.var = 'r2_pre_wm'))

cat('\n=== D2: same table POST-repair ===\n')
print(dcast(tab_post, candidate ~ w_eur, value.var = 'r2_post_wm'))

cat('\n=== D2: pre-repair min eigenvalue (median across blocks) by candidate × w ===\n')
me_tab <- d2[, .(min_eig_pre_med  = median(min_eig_pre,  na.rm = TRUE),
                 min_eig_post_med = median(min_eig_post, na.rm = TRUE)),
             by = .(candidate, w_eur)]
print(dcast(me_tab, candidate ~ w_eur, value.var = 'min_eig_pre_med'))

cat('\n=== D2: rank correlation between (block-mean) pre-repair min-eig and pre-repair M2 across candidates, per w ===\n')
d2_bm <- d2[, .(min_eig_pre = median(min_eig_pre), r2_pre = median(r2_pre)),
            by = .(candidate, w_eur)]
for (w in w_grid) {
  x <- d2_bm[w_eur == w]
  rho <- suppressWarnings(cor(x$min_eig_pre, x$r2_pre, method = 'spearman'))
  cat(sprintf('  w=%.2f: Spearman(min_eig_pre, r2_pre) across %d candidates = %+.3f\n',
              w, nrow(x), rho))
}

cat('\n=== D3: placebo M2 vs MELD-λ1 M2 (POST-repair, δ=0.01, size-weighted mean) by w ===\n')
d3_01 <- d3[delta == 0.01,
            .(placebo_r2_mean = mean(r2, na.rm = TRUE)),
            by = .(w_eur, perm, chr, block_id, m)]
# aggregate per w
d3_agg <- d3_01[, .(placebo_r2_wm = sw(placebo_r2_mean, m)), by = .(w_eur, perm)]
d3_stats <- d3_agg[, .(placebo_mean = mean(placebo_r2_wm),
                        placebo_sd   = sd(placebo_r2_wm),
                        n_perms      = .N), by = w_eur]
l1_stats <- d1_01[candidate == 'MELD_lambda1',
                  .(lambda1_r2_wm_post = sw(r2_post, m)), by = w_eur]
placebo_tab <- merge(d3_stats, l1_stats, by = 'w_eur')
print(placebo_tab)

cat('\n=== D4: λ sweep, POST-repair, δ=0.01, size-weighted mean r² by λ × w ===\n')
d4_01 <- d4[delta == 0.01, .(r2_wm = sw(r2, m)), by = .(w_eur, lambda)]
print(dcast(d4_01, lambda ~ w_eur, value.var = 'r2_wm'))

cat('\n=== D5: (M2_λ1 − M2_λ0) by |ΔAF| tertile × w, size-weighted mean over blocks ===\n')
d5_01 <- d5[delta == 0.01, .(diff_wm = sw(diff, m_bin)), by = .(w_eur, af_bin)]
print(dcast(d5_01, af_bin ~ w_eur, value.var = 'diff_wm'))

cat('\nDONE\n')
