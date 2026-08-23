#!/usr/bin/env Rscript
# gbmi_r10_composition.R — composition scalars per sweep point.
#
# For each (trait, sweep_id) in gbmi_r10_sweep_points.csv, compute:
#   eff_groups  = 1 / Σ w_p²
#   mean_B_diag = mean over shared chr22 variants of 4·Σ_p w_p (f_p,i − f̄_i)²
# f_p,i comes from the harmonised RDS's per-arm af_{POP} column.
#
# Emits gbmi_r10_composition.csv.

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages(library(data.table))
MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r10')

sw <- fread(file.path(OUT_DIR, 'gbmi_r10_sweep_points.csv'))
TRAITS <- unique(sw$trait)

rows <- list()
for (TR in TRAITS) {
  d <- readRDS(r8_harmonised(TR))
  sw_tr <- sw[trait == TR]

  # Wide weights: sweep_id x pop
  W <- dcast(sw_tr, sweep_kind + sweep_id ~ pop, value.var = 'w', fill = 0)

  # AF matrix (n_snps x n_pops) — arm pops only
  arm_pops <- setdiff(colnames(W), c('sweep_kind','sweep_id'))
  af_cols <- paste0('af_', arm_pops)
  present_af <- af_cols %in% names(d)
  arm_pops_here <- arm_pops[present_af]
  af_cols_here <- af_cols[present_af]
  F_mat <- as.matrix(d[, ..af_cols_here])
  colnames(F_mat) <- arm_pops_here

  for (r in seq_len(nrow(W))) {
    w_vec <- as.numeric(W[r, ..arm_pops_here])
    names(w_vec) <- arm_pops_here
    w_pos <- w_vec[w_vec > 0]
    eff_groups <- 1 / sum(w_pos^2)

    # mean_B_diag needs at least 2 non-zero pops
    if (length(w_pos) >= 2L) {
      Fw <- F_mat[, names(w_pos), drop = FALSE]
      w_norm <- w_pos / sum(w_pos)
      fbar <- as.vector(Fw %*% w_norm)
      dF <- Fw - fbar
      B_diag <- 4 * as.vector((dF^2) %*% w_norm)
      mean_B <- mean(B_diag, na.rm = TRUE)
    } else {
      mean_B <- 0
    }
    rows[[length(rows) + 1L]] <- data.table(
      trait = TR,
      sweep_kind = W$sweep_kind[r],
      sweep_id   = W$sweep_id[r],
      eff_groups = eff_groups,
      mean_B_diag = mean_B,
      n_pops_positive = length(w_pos),
      top_pop_weight = max(w_vec)
    )
  }
}

cs <- rbindlist(rows)
fwrite(cs, file.path(OUT_DIR, 'gbmi_r10_composition.csv'))
cat(sprintf('wrote %s (%d rows)\n',
            file.path(OUT_DIR, 'gbmi_r10_composition.csv'), nrow(cs)))
cat('eff_groups range:\n')
print(range(cs$eff_groups))
