#!/usr/bin/env Rscript
# eval_meld_ld_r10.R — Round 10 Section 2/3 sweep runner.
#
# For one trait, iterate every sweep point in gbmi_r10_sweep_points.csv:
#   - Construct z_synth per SNP from arm β/se rescaled to N_tilde = sweep weights
#     (multiplied by a large constant so absolute magnitudes are stable; the
#     construction is scale-invariant in N_tilde so the multiplier drops out).
#   - For each block: sub-index to variants present in the harmonised RDS,
#     build candidate LD panels using reconstruct_R_lambda_P at sweep weights,
#     equal weights, λ=1 at sweep weights, each single-pop, and one scrambled-
#     weight vector. Apply psd_repair to each. Run M2 (5 mask draws × 3 δ).
#   - Block bootstrap CI per candidate × δ; paired diffs against MELD-λ0-true.
#
# Emits gbmi_r10_m2_{per_block,ci,paired_ci,eigen}.csv with sweep_kind and
# sweep_id columns. Appends by trait so a partial rerun replaces only that
# trait's rows.
#
# Usage: Rscript eval_meld_ld_r10.R <trait>

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
MELD_LD_DIR <- file.path(MISC, 'meld_ld')
source(file.path(MISC, 'm2_core.R'))
source(file.path(MISC, 'gbmi_r10_synth.R'))

args  <- commandArgs(trailingOnly = TRUE)
TRAIT <- args[1]
if (is.na(TRAIT)) stop('usage: eval_meld_ld_r10.R <trait>')

CHR          <- 22L
POPS         <- c('EUR','EAS','AFR','CSA','AMR')
MASK_FRAC    <- 0.10
N_DRAWS_MASK <- 3L
DELTAS       <- c(0.001, 0.01, 0.1)
B_BOOT       <- 2000L
SEED         <- 10000L
MIN_BLOCK_M  <- 30L

# Load sweep points, filter to this trait
sw_all_pts <- fread(file.path(MISC, 'gbmi_r10_sweep_points.csv'))
pts <- sw_all_pts[trait == TRAIT]
if (!nrow(pts)) stop(sprintf('no sweep points for %s', TRAIT))
sweep_ids <- unique(pts[, .(sweep_kind, sweep_id)])
cat(sprintf('=== R10 sweep: %s (%d sweep points) ===\n', TRAIT, nrow(sweep_ids)))

# Harmonised trait RDS
d_all <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', TRAIT)))
# Which pops have arms for this trait?
arm_pops <- intersect(POPS, unique(pts$pop))
cat('arm pops:', paste(arm_pops, collapse=', '), '\n')

# Load block RDSes once
chr_dir <- file.path(MELD_LD_DIR, sprintf('chr%d', CHR))
block_files <- sort(list.files(chr_dir, pattern = '^block_.*\\.rds$', full.names = TRUE))
cat(sprintf('chr%d: %d blocks (loading into memory)\n', CHR, length(block_files)))
blocks <- lapply(block_files, readRDS)

# Precompute per-block index into d_all
block_slice <- list()
for (bi in seq_along(blocks)) {
  bl <- blocks[[bi]]
  in_ss <- bl$SNP %in% d_all$rsid
  if (sum(in_ss) < MIN_BLOCK_M) next
  idx <- which(in_ss)
  slice <- list(SNP = bl$SNP[idx])
  for (P in POPS) {
    slice[[paste0('R_', P)]] <- bl[[paste0('R_', P)]][idx, idx, drop = FALSE]
    slice[[paste0('v_', P)]] <- bl[[paste0('v_', P)]][idx]
    slice[[paste0('f_', P)]] <- bl[[paste0('f_', P)]][idx]
  }
  block_slice[[as.character(bl$block_id)]] <- list(
    block_id = bl$block_id, slice = slice, m = length(idx))
}
cat(sprintf('kept %d/%d blocks with M ≥ %d\n', length(block_slice), length(blocks), MIN_BLOCK_M))

# Precompute per-block d_all rows (rsid → row index)
d_rowmap <- setNames(seq_len(nrow(d_all)), d_all$rsid)

# Function: build weight vector over POPS from a sweep_id row set
pts_wide <- dcast(pts, sweep_kind + sweep_id ~ pop, value.var = 'w', fill = 0)
pts_wide[, one_line_id := paste0(sweep_kind, ':', sweep_id)]

# Deterministic scrambled-weight permutation per sweep point (seeded by hash)
scramble <- function(w_vec, salt = 0L) {
  set.seed(as.integer(sum(utf8ToInt(paste0(w_vec, collapse = '_'))) + salt))
  ord <- sample(length(w_vec))
  setNames(w_vec[ord], names(w_vec))
}

m2_rows  <- list()
eig_rows <- list()

for (si in seq_len(nrow(pts_wide))) {
  row <- pts_wide[si]
  kind <- row$sweep_kind
  sid  <- row$sweep_id
  # Weight vector over POPS (may include zero for non-arm pops in this trait)
  w_true <- setNames(rep(0, length(POPS)), POPS)
  for (P in POPS) {
    v <- row[[P]]
    if (length(v) && is.finite(v)) w_true[P] <- v
  }
  # Sanity: sum to 1
  stopifnot(abs(sum(w_true) - 1) < 1e-9)

  # Sweep pops actually contributing
  w_arm_pos <- w_true[w_true > 0 & names(w_true) %in% arm_pops]
  if (length(w_arm_pos) < 1L) next

  # Construct synth z-scores. N_tilde vector uses a fixed scale (10000 * w) —
  # scale is irrelevant since z is scale-invariant, but keeps ratios sane.
  N_tilde <- setNames(rep(0, length(POPS)), POPS)
  N_tilde[names(w_arm_pos)] <- 10000 * w_arm_pos
  N_tilde <- N_tilde[N_tilde > 0]  # drop zero pops for r10_synth
  synth <- r10_synth(d_all, N_tilde)
  # SNP → z map
  z_by_snp <- data.table(SNP = synth$rsid, z = synth$z_synth)

  # Equal-weights vector (over arm pops of this trait only)
  w_equal <- setNames(rep(0, length(POPS)), POPS)
  w_equal[arm_pops] <- 1 / length(arm_pops)

  # Scrambled — permute the sweep weights across arm pops
  w_scr <- setNames(rep(0, length(POPS)), POPS)
  w_arm_vec <- setNames(w_true[arm_pops], arm_pops)
  w_scr[arm_pops] <- as.numeric(scramble(w_arm_vec, salt = si))

  # Which single-pop candidates? All arm pops present in the block set.
  single_pops <- arm_pops

  # Iterate blocks
  for (bi in seq_along(block_slice)) {
    bs <- block_slice[[bi]]
    slice <- bs$slice
    m <- bs$m
    z_local <- z_by_snp[match(slice$SNP, SNP), z]
    if (any(!is.finite(z_local))) next

    # Build candidates
    cand_pre <- list()
    cand_pre[['MELD_lambda0_true']]  <- reconstruct_R_lambda_P(slice, 0, w_true,  POPS)
    cand_pre[['MELD_lambda0_equal']] <- reconstruct_R_lambda_P(slice, 0, w_equal, POPS)
    cand_pre[['MELD_lambda1_true']]  <- reconstruct_R_lambda_P(slice, 1, w_true,  POPS)
    cand_pre[['MELD_lambda0_scrambled']] <- reconstruct_R_lambda_P(slice, 0, w_scr,  POPS)
    for (P in single_pops) cand_pre[[P]] <- slice[[paste0('R_', P)]]

    cand_R <- lapply(cand_pre, psd_repair)

    for (cn in names(cand_R)) {
      eig_rows[[length(eig_rows) + 1L]] <- data.table(
        trait = TRAIT, sweep_kind = kind, sweep_id = sid,
        chr = CHR, block_id = bs$block_id, m = m,
        candidate = cn,
        min_eig_pre  = min_eig(cand_pre[[cn]]),
        min_eig_post = min_eig(cand_R[[cn]])
      )
    }

    set.seed(SEED + si * 137L + bi)
    for (dr in seq_len(N_DRAWS_MASK)) {
      mask_S <- sort(sample.int(m, size = max(2L, round(MASK_FRAC * m))))
      for (delta in DELTAS) {
        for (cn in names(cand_R)) {
          r2 <- run_M2_one(cand_R[[cn]], z_local, mask_S, delta)
          m2_rows[[length(m2_rows) + 1L]] <- data.table(
            trait = TRAIT, sweep_kind = kind, sweep_id = sid,
            chr = CHR, block_id = bs$block_id, m = m,
            candidate = cn, delta = delta, mask_draw = dr, r2 = r2
          )
        }
      }
    }
  }
  if (si %% 5L == 0L) cat(sprintf('  %d/%d sweep points done\n', si, nrow(pts_wide)))
}
cat(sprintf('done %d sweep points\n', nrow(pts_wide)))

m2  <- rbindlist(m2_rows)
eig <- rbindlist(eig_rows)

# Average mask draws within block
m2_block <- m2[, .(r2 = mean(r2, na.rm = TRUE), m = first(m)),
               by = .(trait, sweep_kind, sweep_id, chr, block_id, candidate, delta)]

# Block bootstrap CI per (sweep_id, candidate, delta)
set.seed(SEED + 999L)
ci_rows <- list()
for (sid in unique(m2_block$sweep_id)) {
  for (dl in DELTAS) {
    for (cn in unique(m2_block$candidate)) {
      sub <- m2_block[sweep_id == sid & candidate == cn & delta == dl]
      if (nrow(sub) == 0L) next
      st <- boot_mean(sub$r2, sub$m, B = B_BOOT)
      ci_rows[[length(ci_rows) + 1L]] <- data.table(
        trait = TRAIT,
        sweep_kind = sub$sweep_kind[1], sweep_id = sid,
        delta = dl, candidate = cn,
        mean = st[['mean']], sd = st[['sd']],
        lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
      )
    }
  }
}
ci <- rbindlist(ci_rows)

# Paired diffs against MELD-λ0-true, plus λ1-λ0 and λ0-true-λ0-equal (basin)
paired_rows <- list()
m2_wide <- dcast(m2_block,
                 trait + sweep_kind + sweep_id + chr + block_id + m + delta ~ candidate,
                 value.var = 'r2')
for (sid in unique(m2_wide$sweep_id)) {
  for (dl in DELTAS) {
    sub <- m2_wide[sweep_id == sid & delta == dl]
    base <- 'MELD_lambda0_true'
    comparisons <- list()
    for (P in intersect(arm_pops, names(sub))) comparisons[[length(comparisons)+1]] <- c(base, P)
    comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_equal')       # §4.3
    comparisons[[length(comparisons)+1]] <- c(base, 'MELD_lambda0_scrambled')    # scrambled check
    comparisons[[length(comparisons)+1]] <- c('MELD_lambda1_true', base)         # §4.2
    for (pn in comparisons) {
      if (!(pn[1] %in% names(sub) && pn[2] %in% names(sub))) next
      a <- sub[[pn[1]]]; b <- sub[[pn[2]]]
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 5L) next
      st <- boot_paired(a[ok], b[ok], sub$m[ok], B = B_BOOT)
      paired_rows[[length(paired_rows) + 1L]] <- data.table(
        trait = TRAIT,
        sweep_kind = sub$sweep_kind[1], sweep_id = sid,
        delta = dl,
        comparison = sprintf('%s − %s', pn[1], pn[2]),
        mean_diff = st[['mean_diff']], sd = st[['sd']],
        lo95 = st[['lo95.2.5%']], hi95 = st[['hi95.97.5%']]
      )
    }
  }
}
paired <- rbindlist(paired_rows)

append_or_start <- function(dt, fp) {
  if (file.exists(fp)) {
    prior <- fread(fp)
    prior <- prior[trait != TRAIT]
    dt <- rbind(prior, dt, use.names = TRUE, fill = TRUE)
  }
  fwrite(dt, fp)
}
append_or_start(m2_block, file.path(MISC, 'gbmi_r10_m2_per_block.csv'))
append_or_start(ci,       file.path(MISC, 'gbmi_r10_m2_ci.csv'))
append_or_start(paired,   file.path(MISC, 'gbmi_r10_m2_paired_ci.csv'))
append_or_start(eig,      file.path(MISC, 'gbmi_r10_m2_eigen.csv'))

cat(sprintf('\nwrote per_block (%d) ci (%d) paired (%d) eigen (%d) rows\n',
            nrow(m2_block), nrow(ci), nrow(paired), nrow(eig)))
cat('DONE\n')
