#!/usr/bin/env Rscript
# gbmi_r10_sweep_config.R — enumerate sweep points per trait.
#
# Emits long-form gbmi_r10_sweep_points.csv:
#   trait, sweep_kind, sweep_id, pop, w
# with w summing to 1 across pops within (trait, sweep_kind, sweep_id).
#
# sweep_kind ∈ {'EUR_EAS', 'EUR_AFR', 'multi'}:
#   EUR_EAS: 11 points, w_EUR ∈ {0, 0.1, ..., 1.0}
#   EUR_AFR: 11 points, w_EUR ∈ {0, 0.1, ..., 1.0}
#   multi:   4 named points using all arms in the trait:
#     'equal'   — 1/K each
#     'reported'— proportional to reported per-arm median N
#     'w50rest' — EUR = 0.5, remainder split evenly across non-EUR arms
#     'w25rest' — EUR = 0.25, remainder split evenly across non-EUR arms

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages(library(data.table))
MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'gbmi_r10_synth.R'))

TRAITS <- c('Asthma','COPD','Gout','HF','IPF','Stroke','VTE')
ARM_POPS <- list(
  Asthma = c('EUR','EAS','AFR','CSA','AMR'),
  COPD   = c('EUR','EAS','AFR','AMR'),
  Gout   = c('EUR','EAS','AFR','AMR'),
  HF     = c('EUR','EAS','AFR','AMR'),
  IPF    = c('EUR','EAS','AFR','AMR'),
  Stroke = c('EUR','EAS','AFR','AMR'),
  VTE    = c('EUR','EAS','AFR','AMR')
)

W_EUR_GRID <- seq(0, 1, by = 0.1)

rows <- list()
push <- function(trait, kind, id, pops, w) {
  # Normalise (defensive)
  w <- w / sum(w)
  for (i in seq_along(pops)) {
    rows[[length(rows) + 1L]] <<- data.table(
      trait = trait, sweep_kind = kind, sweep_id = id,
      pop = pops[i], w = w[i]
    )
  }
}

for (TR in TRAITS) {
  pops <- ARM_POPS[[TR]]
  d <- readRDS(file.path(MISC, sprintf('gbmi_r8_%s_chr22.rds', TR)))

  # --- EUR/EAS sweep ---
  for (k in seq_along(W_EUR_GRID)) {
    w_eur <- W_EUR_GRID[k]
    w <- setNames(rep(0, length(pops)), pops)
    w['EUR'] <- w_eur
    w['EAS'] <- 1 - w_eur
    push(TR, 'EUR_EAS', sprintf('EUR_EAS_%02d', k), pops, w)
  }

  # --- EUR/AFR sweep ---
  for (k in seq_along(W_EUR_GRID)) {
    w_eur <- W_EUR_GRID[k]
    w <- setNames(rep(0, length(pops)), pops)
    w['EUR'] <- w_eur
    w['AFR'] <- 1 - w_eur
    push(TR, 'EUR_AFR', sprintf('EUR_AFR_%02d', k), pops, w)
  }

  # --- Multi-pop points ---
  K <- length(pops)
  # equal
  w <- setNames(rep(1 / K, K), pops)
  push(TR, 'multi', 'equal', pops, w)
  # reported
  N_rep <- r10_reported_N(d, pops)
  w_rep <- N_rep / sum(N_rep)
  push(TR, 'multi', 'reported', pops, w_rep)
  # w50rest — EUR=0.5, rest evenly across non-EUR
  non_eur <- setdiff(pops, 'EUR')
  w <- setNames(rep(0, K), pops)
  w['EUR'] <- 0.5
  w[non_eur] <- 0.5 / length(non_eur)
  push(TR, 'multi', 'w50rest', pops, w)
  # w25rest — EUR=0.25, rest evenly
  w <- setNames(rep(0, K), pops)
  w['EUR'] <- 0.25
  w[non_eur] <- 0.75 / length(non_eur)
  push(TR, 'multi', 'w25rest', pops, w)
}

sw <- rbindlist(rows)
fwrite(sw, file.path(MISC, 'gbmi_r10_sweep_points.csv'))
cat(sprintf('wrote %s (%d rows)\n',
            file.path(MISC, 'gbmi_r10_sweep_points.csv'), nrow(sw)))
# Sanity: sum of w per (trait, sweep_id) should be 1
chk <- sw[, .(sumw = sum(w)), by = .(trait, sweep_kind, sweep_id)]
stopifnot(all(abs(chk$sumw - 1) < 1e-9))
cat('per-point weights sum to 1 ✓\n')
cat(sprintf('  points per trait: EUR_EAS=%d, EUR_AFR=%d, multi=%d\n',
            length(W_EUR_GRID), length(W_EUR_GRID), 4L))
cat(sprintf('  total points across 7 traits: %d\n',
            nrow(unique(sw[, .(trait, sweep_id)]))))
