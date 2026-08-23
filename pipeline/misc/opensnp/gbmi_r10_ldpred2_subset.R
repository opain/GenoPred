#!/usr/bin/env Rscript
# gbmi_r10_ldpred2_subset.R — Round 10 Section 5.2 LDpred2-auto subset.
#
# For a chosen subset of (trait, sweep_id, candidate), invoke
# bigsnpr::snp_ldpred2_auto directly on:
#   - constructed IVW meta sumstats at the sweep weights (β_synth, se_synth,
#     nominal N_eff = sum of N_tilde)
#   - candidate block-diagonal R matrix (MELD-λ₀-true, MELD-λ₀-equal, EUR,
#     or EAS) restricted to chr22 and matched to the sumstats.
#
# Records: n retained chains, mean h²_est across retained chains,
# mean p_est, log10-scale ratio, convergence flag.
#
# Scope (per plan):
#   2 traits × 3 EUR/EAS sweep points × 4 candidates = 24 runs.
#
# Bypasses the pipeline's ldpred2.R since we don't need PLINK genotypes or
# a full refdir — just the block-diagonal correlation matrix in memory.

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(Matrix)
  library(bigsparser)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
MELD_LD_DIR <- file.path(MISC, 'meld_ld', 'chr22')
source(file.path(MISC, 'm2_core.R'))
source(file.path(MISC, 'gbmi_r10_synth.R'))

POPS_ALL <- c('EUR','EAS','AFR','CSA','AMR')

# Subset scope
TRAITS_SUB <- c('Asthma', 'COPD')
# EUR_EAS_01 = w_EUR=0, EUR_EAS_06 = w_EUR=0.5, EUR_EAS_11 = w_EUR=1.0
SWEEP_IDS_SUB <- c('EUR_EAS_01', 'EUR_EAS_06', 'EUR_EAS_11')
CANDIDATES <- c('MELD_lambda0_true', 'MELD_lambda0_equal', 'EUR', 'EAS')
SEED <- 20000L
N_ITER <- 500L
N_BURN <- 200L
N_CHAIN <- 20L

sw_all <- fread(file.path(MISC, 'gbmi_r10_sweep_points.csv'))

# Load R8 blocks once
block_files <- sort(list.files(MELD_LD_DIR, pattern = '^block_.*\\.rds$', full.names = TRUE))
blocks <- lapply(block_files, readRDS)

build_candidate_block_R <- function(bl, cand, w_true, arm_pops) {
  # Build a full-block R (all its SNPs, no sub-index) for a given candidate.
  slice <- list(SNP = bl$SNP)
  for (P in POPS_ALL) {
    slice[[paste0('R_', P)]] <- bl[[paste0('R_', P)]]
    slice[[paste0('v_', P)]] <- bl[[paste0('v_', P)]]
    slice[[paste0('f_', P)]] <- bl[[paste0('f_', P)]]
  }
  R <- if (cand == 'MELD_lambda0_true') {
    reconstruct_R_lambda_P(slice, 0, w_true, POPS_ALL)
  } else if (cand == 'MELD_lambda0_equal') {
    w_eq <- setNames(rep(0, length(POPS_ALL)), POPS_ALL)
    w_eq[arm_pops] <- 1 / length(arm_pops)
    reconstruct_R_lambda_P(slice, 0, w_eq, POPS_ALL)
  } else if (cand %in% POPS_ALL) {
    slice[[paste0('R_', cand)]]
  } else stop('unknown candidate ', cand)
  psd_repair(R)
}

# Assemble the block-diagonal sparse R and per-block metadata for the chr,
# restricted to variants present in the sumstats data (by rsid).
build_chr_corr <- function(cand_name, w_true, arm_pops, sumstats_rsids) {
  block_R_list <- list()
  ids_ord <- character(0)
  cum <- 0L
  block_intervals <- data.table(block_id = integer(0),
                                start = integer(0),
                                end = integer(0))
  for (bl in blocks) {
    common_ix <- which(bl$SNP %in% sumstats_rsids)
    if (length(common_ix) < 2L) next
    slice <- list(SNP = bl$SNP[common_ix])
    for (P in POPS_ALL) {
      slice[[paste0('R_', P)]] <- bl[[paste0('R_', P)]][common_ix, common_ix, drop = FALSE]
      slice[[paste0('v_', P)]] <- bl[[paste0('v_', P)]][common_ix]
      slice[[paste0('f_', P)]] <- bl[[paste0('f_', P)]][common_ix]
    }
    R <- if (cand_name == 'MELD_lambda0_true') {
      reconstruct_R_lambda_P(slice, 0, w_true, POPS_ALL)
    } else if (cand_name == 'MELD_lambda0_equal') {
      w_eq <- setNames(rep(0, length(POPS_ALL)), POPS_ALL)
      w_eq[arm_pops] <- 1 / length(arm_pops)
      reconstruct_R_lambda_P(slice, 0, w_eq, POPS_ALL)
    } else if (cand_name %in% POPS_ALL) {
      slice[[paste0('R_', cand_name)]]
    } else stop('unknown candidate ', cand_name)
    R <- psd_repair(R)
    block_R_list[[length(block_R_list) + 1L]] <- as(R, 'sparseMatrix')
    m <- nrow(R)
    block_intervals <- rbind(block_intervals,
                             data.table(block_id = bl$block_id,
                                        start = cum + 1L,
                                        end = cum + m))
    cum <- cum + m
    ids_ord <- c(ids_ord, bl$SNP[common_ix])
  }
  corr_bd <- as(Matrix::bdiag(block_R_list), 'symmetricMatrix')
  attr(corr_bd, 'rsids') <- ids_ord
  attr(corr_bd, 'block_intervals') <- block_intervals
  corr_bd
}

# ---------------------------------------------------------------------------
rows <- list()
for (TR in TRAITS_SUB) {
  d <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', TR)))
  arm_pops <- c('EUR','EAS','AFR','CSA','AMR')
  arm_pops <- arm_pops[arm_pops %in% sub('^N_', '', grep('^N_', names(d), value = TRUE))]

  for (SID in SWEEP_IDS_SUB) {
    # Sweep weights over arm pops
    sp <- sw_all[trait == TR & sweep_id == SID]
    if (!nrow(sp)) next
    w_true <- setNames(rep(0, length(POPS_ALL)), POPS_ALL)
    for (r in seq_len(nrow(sp))) w_true[sp$pop[r]] <- sp$w[r]

    # Constructed sumstats: N_tilde = w × real total meta N per trait.
    # z_synth is scale-invariant; β_synth/se_synth scale with this so
    # LDpred2's n_eff reflects the actual GWAS precision.
    N_real <- r10_reported_N(d, arm_pops)
    N_real_total <- sum(N_real, na.rm = TRUE)
    N_tilde <- N_real_total * w_true[w_true > 0 & names(w_true) %in% arm_pops]
    N_tilde <- N_tilde[N_tilde > 0]
    synth <- r10_synth(d, N_tilde)
    ss <- data.table(rsid = d$rsid,
                     chr  = 22,
                     pos  = d$pos,
                     a0   = d$ref,
                     a1   = d$alt,
                     beta_synth = synth$beta_synth,
                     se_synth   = synth$se_synth,
                     z_synth    = synth$z_synth,
                     ok         = synth$ok)
    ss <- ss[ok == TRUE & is.finite(beta_synth) & is.finite(se_synth) & se_synth > 0]
    # n_eff = sum of N_tilde (nominal; LDpred2 uses it for its posterior calc)
    n_eff <- sum(N_tilde)
    ss[, n_eff := n_eff]

    for (CAND in CANDIDATES) {
      cat(sprintf('\n=== %s | %s | %s ===\n', TR, SID, CAND))
      corr <- build_chr_corr(CAND, w_true, arm_pops, ss$rsid)
      rsids_R <- attr(corr, 'rsids')
      keep <- ss[rsid %in% rsids_R]
      keep <- keep[match(rsids_R, rsid)]
      stopifnot(all(keep$rsid == rsids_R))

      # Initial h² from LDSC-lite: h² = mean(chi²) - 1 divided by mean ld ...
      # Simpler: use variance of z scaled by 1/n_eff. Or seed at h²=0.2.
      h2_init <- max(0.01, min(0.5,
                               (mean(keep$z_synth^2, na.rm = TRUE) - 1) /
                                 (mean(keep$z_synth^2, na.rm = TRUE))))
      cat(sprintf('  n_snps=%d  n_eff=%.0f  h2_init=%.3f\n',
                  nrow(keep), n_eff, h2_init))

      df_beta <- data.frame(beta = keep$beta_synth,
                            beta_se = keep$se_synth,
                            n_eff = keep$n_eff)
      # bigsnpr wants an SFBM, not a dsCMatrix
      corr_sfbm <- tryCatch(
        as_SFBM(as(corr, 'dsCMatrix'), backingfile = tempfile(fileext = '.sfbm')),
        error = function(e) { cat(sprintf('  SFBM ERROR: %s\n', conditionMessage(e))); NULL }
      )
      auto <- if (is.null(corr_sfbm)) NULL else tryCatch({
        snp_ldpred2_auto(
          corr_sfbm, df_beta, h2_init,
          vec_p_init = seq_log(1e-4, 0.5, N_CHAIN),
          burn_in = N_BURN, num_iter = N_ITER,
          ncores = 4L, sparse = FALSE
        )
      }, error = function(e) {
        cat(sprintf('  ERROR: %s\n', conditionMessage(e)))
        NULL
      })

      if (is.null(auto)) {
        rows[[length(rows) + 1L]] <- data.table(
          trait = TR, sweep_id = SID, candidate = CAND,
          n_chains = N_CHAIN, chains_kept = 0L,
          h2_est = NA_real_, h2_sd = NA_real_,
          p_est = NA_real_, p_sd = NA_real_,
          converged = FALSE, n_snps = nrow(keep), n_eff = n_eff
        )
        next
      }

      # Chain-keep filter per meld_context_2 / ldpred2.R:297-299
      ranges <- sapply(auto, function(a) diff(range(a$corr_est)))
      threshold <- 0.95 * quantile(ranges, 0.95, na.rm = TRUE)
      keep_ch <- which(ranges >= threshold)

      h2_vec <- sapply(auto, `[[`, 'h2_est')
      p_vec  <- sapply(auto, `[[`, 'p_est')

      rows[[length(rows) + 1L]] <- data.table(
        trait = TR, sweep_id = SID, candidate = CAND,
        n_chains = length(auto), chains_kept = length(keep_ch),
        h2_est = mean(h2_vec[keep_ch], na.rm = TRUE),
        h2_sd  = sd  (h2_vec[keep_ch], na.rm = TRUE),
        p_est  = mean(p_vec[keep_ch], na.rm = TRUE),
        p_sd   = sd  (p_vec[keep_ch], na.rm = TRUE),
        converged = length(keep_ch) >= 2L,
        n_snps = nrow(keep), n_eff = n_eff
      )
      cat(sprintf('  chains_kept=%d/%d  h2=%.3f  p=%.4f\n',
                  length(keep_ch), length(auto),
                  mean(h2_vec[keep_ch], na.rm = TRUE),
                  mean(p_vec[keep_ch], na.rm = TRUE)))
    }
  }
}

out <- rbindlist(rows)
fwrite(out, file.path(MISC, 'gbmi_r10_ldpred2_subset.csv'))
cat(sprintf('\nwrote %s (%d rows)\n',
            file.path(MISC, 'gbmi_r10_ldpred2_subset.csv'), nrow(out)))
