# gbmi_r10_synth.R — Round 10 Section 1 constructor.
#
# For a given trait's harmonised RDS (from Round 8) and a nominal per-arm
# sample-size vector `N_tilde`, construct per-SNP β/se/z of the IVW meta
# under the chosen composition. The construction uses the arms' actual
# β and se — not `2f(1-f)` reference-panel proxies — so the SEs carry
# real per-SNP precision, imputation quality and case/control structure.
#
# Formula (per prompt §1):
#   w_tilde_p,i = (N_tilde_p / N_p) * (1 / se_p,i^2)
#   z_synth,i   = sum_p sqrt(w_tilde_p,i) * z_p,i / sqrt(sum_p w_tilde_p,i)
#   beta_synth,i = sum_p w_tilde_p,i * beta_p,i / sum_p w_tilde_p,i
#   se_synth,i   = 1 / sqrt(sum_p w_tilde_p,i)
#
# N_p is the arm's median per-SNP N.
#
# Sourced by gbmi_r10_sanity.R and eval_meld_ld_r10.R. Pure functions, no
# side effects on load.

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Given a harmonised trait RDS and a named N_tilde vector (over arm pops
# present in `d`), return a data.table with rsid + beta_synth + se_synth
# + z_synth. Drops rows where any pop's se is NA/0/inf.
#
# `d` must carry columns beta_{POP}, se_{POP}, N_{POP} for every POP name
# in `N_tilde`. Ignored pops (N_tilde == 0) don't need columns.

r10_synth <- function(d, N_tilde) {
  stopifnot(is.data.table(d))
  pops <- names(N_tilde)[N_tilde > 0]
  stopifnot(length(pops) >= 1L)

  # Sanity: required arm columns present
  need <- c(sprintf('beta_%s', pops),
            sprintf('se_%s',   pops),
            sprintf('N_%s',    pops))
  missing <- setdiff(need, names(d))
  if (length(missing) > 0)
    stop(sprintf('missing arm columns: %s', paste(missing, collapse = ', ')))

  # Median per-SNP N per arm — the denominator N_p
  N_p <- sapply(pops, function(p) as.numeric(median(d[[sprintf('N_%s', p)]], na.rm = TRUE)))
  names(N_p) <- pops

  # Rescale factor per pop
  ratio <- N_tilde[pops] / N_p[pops]  # scalar per pop

  # Build matrices of β, se, z per (SNP, pop). SNP-major.
  B <- do.call(cbind, lapply(pops, function(p) d[[sprintf('beta_%s', p)]]))
  S <- do.call(cbind, lapply(pops, function(p) d[[sprintf('se_%s',   p)]]))
  colnames(B) <- colnames(S) <- pops
  Z <- B / S

  # Rescaled precision per (SNP, pop): w_tilde_p,i = ratio_p * (1/se_p,i^2)
  W <- t(t(1 / (S^2)) * ratio)  # multiply column p by ratio[p]

  # Zero out entries where SE is not finite/positive (bad arm rows)
  bad <- !is.finite(S) | (S <= 0) | !is.finite(B)
  W[bad] <- 0

  sumW <- rowSums(W)
  # Drop rows with no contributing pop
  ok <- sumW > 0

  beta_synth <- rowSums(W * B, na.rm = TRUE) / sumW
  se_synth   <- 1 / sqrt(sumW)
  z_synth    <- rowSums(sqrt(W) * Z, na.rm = TRUE) / sqrt(sumW)

  out <- data.table(rsid = d$rsid,
                    beta_synth = beta_synth,
                    se_synth   = se_synth,
                    z_synth    = z_synth,
                    ok         = ok)
  # Convenience: N_eff for LDpred2 downstream (Willer 2010 IVW-N-eff)
  # sum of N_tilde is the natural nominal N_eff for this meta.
  attr(out, 'N_tilde') <- N_tilde
  attr(out, 'N_p') <- N_p
  attr(out, 'ratio') <- ratio
  attr(out, 'n_ok') <- sum(ok)
  out
}

# ---------------------------------------------------------------------------
# Convenience: real reported per-arm N vector for a trait, derived from the
# harmonised RDS itself (medians of N_{POP} columns).
r10_reported_N <- function(d, pops) {
  setNames(sapply(pops, function(p) {
    col <- sprintf('N_%s', p)
    if (!(col %in% names(d))) return(NA_real_)
    as.numeric(median(d[[col]], na.rm = TRUE))
  }), pops)
}
