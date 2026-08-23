#!/usr/bin/env Rscript
# noise_inject_m2.R — Round 10 Section 5.1: is M2 sensitive to panel error?
#
# Add calibrated Gaussian noise to a well-estimated EUR panel and re-score
# M2 against a fixed reference trait (Asthma z-scores). If M2 barely moves
# across a 100x range of simulated panel n, then M2 is insensitive to
# panel-estimation error, and every M2 result to date needs re-interpreting.
#
# Baselines used:
#   (a) UKB EUR from meld_ld_ukb/chr22/ (67 sub-blocks, N=362063) —
#       genuinely well-estimated.
#   (b) 1KG+HGDP EUR from meld_ld/chr22/ (24 blocks, N=665) — the R8 panel
#       whose flatness under UKB replacement was one of the motivating
#       observations.
#
# Noise model (per prompt §5.1):
#   For each off-diagonal r_ij: add Gaussian noise with SD = (1 - r_ij^2) / sqrt(n)
#   Then symmetrise (r_ij_noisy + r_ji_noisy)/2, re-apply psd_repair.
#   Diagonal stays at 1.
#
# Nominal n grid: {500, 1000, 5000, 50000}. 5 draws per level.
#
# Emits: gbmi_r10_noise_m2.csv (per block, per baseline, per n, per draw:
#   r2, min_eig_pre, min_eig_post, cond_number).

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

MISC   <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
R8_DIR <- file.path(MISC, 'meld_ld',     'chr22')
UK_DIR <- file.path(MISC, 'meld_ld_ukb', 'chr22')
source(file.path(MISC, 'm2_core.R'))

REF_TRAIT   <- 'Asthma'
POP         <- 'EUR'
NOISE_NS    <- c(500L, 1000L, 5000L, 50000L)
N_DRAWS     <- 5L
MASK_FRAC   <- 0.10
N_DRAWS_MASK <- 5L
DELTAS      <- c(0.001, 0.01, 0.1)
SEED        <- 12000L
MIN_M       <- 30L

d_all <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', REF_TRAIT)))
d_all[, z := meta_beta / meta_se]
z_by_snp <- d_all[, .(SNP = rsid, z)]

perturb_and_repair <- function(R, n) {
  # SD ~ (1 - r^2) / sqrt(n) elementwise, off-diagonals only
  M <- nrow(R)
  if (M < 2L) return(list(R = R, min_eig_pre = 1, min_eig_post = 1, cond = 1))
  E <- matrix(rnorm(M * M) * ((1 - R^2) / sqrt(n)), M, M)
  E[lower.tri(E)] <- t(E)[lower.tri(E)]  # symmetrise
  diag(E) <- 0
  R_noisy <- R + E
  # Clip to [-1, 1] before repair (Gaussian tails can exceed)
  R_noisy[R_noisy > 1] <- 1
  R_noisy[R_noisy < -1] <- -1
  diag(R_noisy) <- 1
  min_pre <- min_eig(R_noisy)
  R_rep <- psd_repair(R_noisy)
  ev_post <- eigen(R_rep, symmetric = TRUE, only.values = TRUE)$values
  min_post <- min(ev_post)
  cond <- if (min(abs(ev_post)) > 0) max(ev_post) / max(min(ev_post), 1e-12) else Inf
  list(R = R_rep, min_eig_pre = min_pre, min_eig_post = min_post, cond = cond)
}

score_M2 <- function(R, z_local, seed_here) {
  set.seed(seed_here)
  m <- length(z_local)
  # Multiple mask draws, mean across draws, mean across deltas at delta=0.01 (headline)
  out <- list()
  for (dr in seq_len(N_DRAWS_MASK)) {
    mask_S <- sort(sample.int(m, size = max(2L, round(MASK_FRAC * m))))
    for (dl in DELTAS) {
      out[[length(out) + 1L]] <- data.table(
        mask_draw = dr, delta = dl,
        r2 = run_M2_one(R, z_local, mask_S, dl))
    }
  }
  rbindlist(out)
}

process_baseline <- function(block_dir, label) {
  block_files <- sort(list.files(block_dir, pattern = '^block_.*\\.rds$', full.names = TRUE))
  cat(sprintf('\n=== baseline = %s (%d blocks) ===\n', label, length(block_files)))
  rows <- list()
  for (fi in seq_along(block_files)) {
    bl <- readRDS(block_files[[fi]])
    R_ref <- bl[[paste0('R_', POP)]]
    if (is.null(R_ref) || nrow(R_ref) < MIN_M) next
    in_ss <- bl$SNP %in% z_by_snp$SNP
    if (sum(in_ss) < MIN_M) next
    idx <- which(in_ss)
    R_sub <- R_ref[idx, idx, drop = FALSE]
    z_local <- z_by_snp[match(bl$SNP[idx], SNP), z]
    if (any(!is.finite(z_local))) next
    m <- length(idx)

    # Baseline (no noise) — but still apply psd_repair for parity
    R_base <- psd_repair(R_sub)
    min_post_base <- min_eig(R_base)
    r2_base <- score_M2(R_base, z_local, SEED + fi)
    for (r in seq_len(nrow(r2_base))) {
      rows[[length(rows) + 1L]] <- data.table(
        baseline = label, block_id = bl$block_id, m = m,
        n_nominal = Inf, draw = 0L,
        mask_draw = r2_base$mask_draw[r], delta = r2_base$delta[r],
        r2 = r2_base$r2[r],
        min_eig_pre = min_eig(R_sub), min_eig_post = min_post_base,
        cond_post = NA_real_
      )
    }

    # Perturbed at each n × draw
    for (n in NOISE_NS) {
      for (draw in seq_len(N_DRAWS)) {
        set.seed(as.integer(SEED + (n %% 100000L) * 7L + draw * 131L + fi))
        pr <- perturb_and_repair(R_sub, n)
        r2_p <- score_M2(pr$R, z_local, SEED + 100L * n + draw + fi)
        for (r in seq_len(nrow(r2_p))) {
          rows[[length(rows) + 1L]] <- data.table(
            baseline = label, block_id = bl$block_id, m = m,
            n_nominal = as.numeric(n), draw = draw,
            mask_draw = r2_p$mask_draw[r], delta = r2_p$delta[r],
            r2 = r2_p$r2[r],
            min_eig_pre = pr$min_eig_pre,
            min_eig_post = pr$min_eig_post,
            cond_post = pr$cond
          )
        }
      }
    }
    if (fi %% 5L == 0L) cat(sprintf('  block %d/%d done\n', fi, length(block_files)))
  }
  rbindlist(rows)
}

all_rows <- rbind(
  process_baseline(R8_DIR, '1KG_HGDP_EUR'),
  process_baseline(UK_DIR, 'UKB_EUR')
)
fwrite(all_rows, file.path(MISC, 'gbmi_r10_noise_m2.csv'))
cat(sprintf('\nwrote %s (%d rows)\n',
            file.path(MISC, 'gbmi_r10_noise_m2.csv'), nrow(all_rows)))

# ---- Summary at delta = 0.01, size-weighted mean across blocks per baseline+n
cat('\n=== headline: mean M2 vs nominal n (delta=0.01, size-weighted, over draws + mask draws) ===\n')
head_dt <- all_rows[delta == 0.01,
                    .(r2_mean = weighted.mean(r2, m, na.rm = TRUE),
                      r2_sd = sd(r2, na.rm = TRUE),
                      min_eig_pre_med = median(min_eig_pre),
                      min_eig_post_med = median(min_eig_post),
                      cond_med = median(cond_post, na.rm = TRUE)),
                    by = .(baseline, n_nominal)]
setorder(head_dt, baseline, n_nominal)
print(head_dt)

cat('\nDONE\n')
