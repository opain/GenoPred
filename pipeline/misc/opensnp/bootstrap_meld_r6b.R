#!/usr/bin/env Rscript
# Round 6b D6: block-bootstrap 95% CIs for the headline candidate/λ differences.
#
# Reads the per-block M2 scores from meld_r6b_d1_m2_all_candidates.csv (which
# holds pre-repair AND post-repair r² at δ=0.001/0.01/0.1) and the λ-sweep
# from meld_r6b_d4_lambda_sweep.csv. Resamples blocks with replacement B times,
# recomputes size-weighted means, and reports bootstrap SD + 95% CIs.
#
# Writes: meld_r6b_d6_bootstrap.csv

suppressPackageStartupMessages(library(data.table))

out_dir <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
d1 <- fread(file.path(out_dir, 'meld_r6b_d1_m2_all_candidates.csv'))
d4 <- fread(file.path(out_dir, 'meld_r6b_d4_lambda_sweep.csv'))
d3 <- fread(file.path(out_dir, 'meld_r6b_d3_placebo.csv'))

B <- 2000L
DELTA <- 0.01
set.seed(2026L)

# ---------------------------------------------------------------------------
# Size-weighted mean helper
sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# Block-bootstrap over a per-block value with a size weight m
boot_mean <- function(vals, weights, n_boot = B) {
  n <- length(vals)
  boots <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    boots[b] <- sw(vals[idx], weights[idx])
  }
  c(mean = sw(vals, weights), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

# Bootstrap for a paired candidate-difference: same block, subtract, then mean.
# Aggregates to per (chr, block_id) row per candidate first.
boot_paired_diff <- function(d_by_block_a, d_by_block_b, m, n_boot = B) {
  # d_by_block_a and _b MUST be aligned by (chr, block_id, w_eur)
  diffs <- d_by_block_a - d_by_block_b
  boots <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(length(diffs), length(diffs), replace = TRUE)
    boots[b] <- sw(diffs[idx], m[idx])
  }
  c(mean_diff = sw(diffs, m), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

# ---------------------------------------------------------------------------
# 1. Marginal r² per (candidate, w, δ, pre/post) with bootstrap CIs
d1_delta <- d1[delta == DELTA]

rows <- list()
for (w in unique(d1_delta$w_eur)) {
  for (cn in unique(d1_delta$candidate)) {
    sub <- d1_delta[w_eur == w & candidate == cn]
    if (nrow(sub) == 0L) next
    pre  <- boot_mean(sub$r2_pre,  sub$m)
    post <- boot_mean(sub$r2_post, sub$m)
    rows[[length(rows) + 1L]] <- data.table(
      w_eur = w, candidate = cn, kind = 'pre',
      mean = pre[['mean']],  sd = pre[['sd']],
      lo95 = pre[['lo95.2.5%']], hi95 = pre[['hi95.97.5%']]
    )
    rows[[length(rows) + 1L]] <- data.table(
      w_eur = w, candidate = cn, kind = 'post',
      mean = post[['mean']],  sd = post[['sd']],
      lo95 = post[['lo95.2.5%']], hi95 = post[['hi95.97.5%']]
    )
  }
}
d6_marginal <- rbindlist(rows)

# ---------------------------------------------------------------------------
# 2. Paired candidate differences with bootstrap CIs (post-repair, δ=0.01)
# Focus: MELD-λ0 vs each single-pop; MELD-λ1 vs MELD-λ0 (residual);
# MELD-λ1 vs placebo (residual).
d1_wide <- dcast(d1_delta, chr + block_id + m + w_eur ~ candidate,
                 value.var = c('r2_pre', 'r2_post'))

pair_rows <- list()
diff_names <- list(
  c('MELD_lambda0', 'AFR'),
  c('MELD_lambda0', 'EUR'),
  c('MELD_lambda1', 'MELD_lambda0'),
  c('MELD_lambda1', 'AFR'),
  c('MELD_lambda1', 'EUR')
)
for (pn in diff_names) {
  a <- paste0('r2_post_', pn[1])
  b <- paste0('r2_post_', pn[2])
  for (w in unique(d1_wide$w_eur)) {
    sub <- d1_wide[w_eur == w & is.finite(get(a)) & is.finite(get(b))]
    if (nrow(sub) == 0L) next
    stats <- boot_paired_diff(sub[[a]], sub[[b]], sub$m)
    pair_rows[[length(pair_rows) + 1L]] <- data.table(
      w_eur = w, comparison = sprintf('%s − %s (post)', pn[1], pn[2]),
      mean_diff = stats[['mean_diff']], sd = stats[['sd']],
      lo95 = stats[['lo95.2.5%']], hi95 = stats[['hi95.97.5%']]
    )
  }
}

# Also add the pre-repair MELD_lambda1 − placebo residual
# Placebo (per perm) → per-block mean over perms → paired diff with MELD_lambda1
d3_bl <- d3[delta == DELTA, .(placebo_r2 = mean(r2, na.rm = TRUE)),
            by = .(chr, block_id, m, w_eur)]
l1_bl <- d1_delta[candidate == 'MELD_lambda1',
                  .(chr, block_id, m, w_eur,
                    lambda1_post = r2_post, lambda1_pre = r2_pre)]
paired <- merge(l1_bl, d3_bl, by = c('chr', 'block_id', 'm', 'w_eur'))
for (w in unique(paired$w_eur)) {
  sub <- paired[w_eur == w & is.finite(lambda1_post) & is.finite(placebo_r2)]
  if (nrow(sub) == 0L) next
  stats <- boot_paired_diff(sub$lambda1_post, sub$placebo_r2, sub$m)
  pair_rows[[length(pair_rows) + 1L]] <- data.table(
    w_eur = w, comparison = 'MELD_lambda1 − placebo (post)',
    mean_diff = stats[['mean_diff']], sd = stats[['sd']],
    lo95 = stats[['lo95.2.5%']], hi95 = stats[['hi95.97.5%']]
  )
}
d6_paired <- rbindlist(pair_rows)

# ---------------------------------------------------------------------------
# 3. λ sweep with bootstrap CIs (post-repair)
lambda_rows <- list()
for (w in unique(d4$w_eur)) {
  for (lam in unique(d4$lambda)) {
    sub <- d4[w_eur == w & lambda == lam & delta == DELTA]
    if (nrow(sub) == 0L) next
    stats <- boot_mean(sub$r2, sub$m)
    lambda_rows[[length(lambda_rows) + 1L]] <- data.table(
      w_eur = w, lambda = lam,
      mean = stats[['mean']],  sd = stats[['sd']],
      lo95 = stats[['lo95.2.5%']], hi95 = stats[['hi95.97.5%']]
    )
  }
}
d6_lambda <- rbindlist(lambda_rows)

fwrite(d6_marginal, file.path(out_dir, 'meld_r6b_d6_marginal_ci.csv'))
fwrite(d6_paired,   file.path(out_dir, 'meld_r6b_d6_paired_ci.csv'))
fwrite(d6_lambda,   file.path(out_dir, 'meld_r6b_d6_lambda_ci.csv'))

# ---------------------------------------------------------------------------
# Console tables
cat('\n=== D6 marginal r² with 95% bootstrap CI (POST-repair, δ=0.01) ===\n')
tab <- d6_marginal[kind == 'post',
                   .(candidate, w_eur,
                     r2 = sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95))]
print(dcast(tab, candidate ~ w_eur, value.var = 'r2'))

cat('\n=== D6 paired difference with 95% bootstrap CI ===\n')
d6_paired[, txt := sprintf('%+.4f [%+.4f, %+.4f]', mean_diff, lo95, hi95)]
print(dcast(d6_paired, comparison ~ w_eur, value.var = 'txt'))

cat('\n=== D6 λ sweep r² with 95% bootstrap CI (post-repair, δ=0.01) ===\n')
d6_lambda[, txt := sprintf('%.3f [%.3f, %.3f]', mean, lo95, hi95)]
print(dcast(d6_lambda, lambda ~ w_eur, value.var = 'txt'))
