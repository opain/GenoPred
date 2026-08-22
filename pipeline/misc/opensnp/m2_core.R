# m2_core.R — Round 8 shared core for M2 evaluation.
#
# Verbatim copies of the Round 7b functions from eval_meld_ld_real.R, with the
# minimum mechanical change needed to sweep δ:
#   - run_M2_one() and boot_mean()/boot_paired() take their tuning constants
#     (delta, B) as arguments rather than reading them from the enclosing scope.
# The maths is unchanged.
#
# Rationale for reuse (from meld_round8_gbmi_prompt.md §5):
#   "unchanged from Round 7b and non-negotiable"
# Source of truth: pipeline/misc/opensnp/eval_meld_ld_real.R:128-276.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# Uniform PSD repair: eigen-clamp + diagonal rescale to unit-diag.
psd_repair <- function(R) {
  ev  <- eigen(R, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  R2  <- ev$vectors %*% (lam * t(ev$vectors))
  d   <- sqrt(pmax(diag(R2), 0)); d[d == 0] <- 1
  R2  <- R2 / outer(d, d); diag(R2) <- 1
  R2
}

# ---------------------------------------------------------------------------
# Analytic P-population MELD target R*_lambda (Wahlund form).
#   Sigma_star = Σ_p w_p R_p * outer(sqrt(v_p), sqrt(v_p))
#   B(i,j)     = 4 Σ_p w_p (f_p,i − f̄_i)(f_p,j − f̄_j),  f̄ = Σ_p w_p f_p
#   Sigma      = Sigma_star + lambda * B
#   R          = Sigma / outer(sqrt(V), sqrt(V))
# `bl` must carry R_{POP}, v_{POP}, f_{POP} for every pop in pops_here.
reconstruct_R_lambda_P <- function(bl, lambda, w_pop, pops_here) {
  m <- length(bl$SNP)
  Sigma_star <- matrix(0, m, m)
  V_star     <- rep(0, m)
  for (p in pops_here) {
    R_p   <- bl[[paste0('R_', p)]]
    v_p   <- bl[[paste0('v_', p)]]
    Cov_p <- R_p * outer(sqrt(v_p), sqrt(v_p))
    Sigma_star <- Sigma_star + w_pop[[p]] * Cov_p
    V_star     <- V_star     + w_pop[[p]] * v_p
  }
  if (abs(lambda) > 1e-12) {
    F   <- do.call(cbind, lapply(pops_here, function(p) bl[[paste0('f_', p)]]))
    fbar <- as.vector(F %*% w_pop[pops_here])
    dF  <- F - fbar
    W   <- diag(w_pop[pops_here])
    B   <- 4 * (dF %*% W %*% t(dF))
    Sigma <- Sigma_star + lambda * B
    V     <- V_star     + lambda * diag(B)
  } else {
    Sigma <- Sigma_star
    V     <- V_star
  }
  inv_sd <- 1 / sqrt(V)
  R <- Sigma * outer(inv_sd, inv_sd)
  R[!is.finite(R)] <- 0
  R
}

# ---------------------------------------------------------------------------
# M2 (masked-z re-imputation) score on a single candidate R × mask.
# Returns r² between predicted and held-out z. Uses ridge δ · mean(|λ_MM|).
run_M2_one <- function(R_cand, z, mask_S, delta) {
  M    <- setdiff(seq_along(z), mask_S)
  R_MM <- R_cand[M, M, drop = FALSE]
  R_SM <- R_cand[mask_S, M, drop = FALSE]
  z_M  <- z[M]; z_S <- z[mask_S]
  ev   <- suppressWarnings(eigen(R_MM, symmetric = TRUE))
  U    <- ev$vectors; lam <- ev$values
  ridge <- delta * abs(sum(lam) / length(lam))
  Utz  <- crossprod(U, z_M)
  z_hat <- drop(R_SM %*% (U %*% (Utz / (lam + ridge))))
  r <- suppressWarnings(cor(z_hat, z_S))
  if (is.finite(r)) r^2 else NA_real_
}

# ---------------------------------------------------------------------------
# Size-weighted mean and paired-difference bootstrap CIs.
sw <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

boot_mean <- function(vals, weights, B) {
  n <- length(vals)
  boots <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    boots[b] <- sw(vals[idx], weights[idx])
  }
  c(mean = sw(vals, weights), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

boot_paired <- function(a, b, w, B) {
  diffs <- a - b
  boots <- numeric(B)
  for (i in seq_len(B)) {
    idx <- sample.int(length(diffs), length(diffs), replace = TRUE)
    boots[i] <- sw(diffs[idx], w[idx])
  }
  c(mean_diff = sw(diffs, w), sd = sd(boots),
    lo95 = quantile(boots, 0.025), hi95 = quantile(boots, 0.975))
}

# ---------------------------------------------------------------------------
# Min-eigenvalue helper — for pre/post-repair reporting per candidate.
min_eig <- function(R) {
  ev <- suppressWarnings(eigen(R, symmetric = TRUE, only.values = TRUE))
  min(ev$values)
}
