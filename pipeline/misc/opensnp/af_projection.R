# af_projection.R — extracted from eval_meld_ld_real.R:63-114.
#
# Ancestry-composition estimator from a GWAS allele-frequency vector, projected
# against Privé's 1KG+HGDP fine-population reference, rolled up to the GenoPred
# super-populations {EUR, EAS, AFR, CSA, AMR, MID}.
#
# Inputs to compute_afproj_weights():
#   gwas_freq: data.table with columns
#     chr (int), pos (int), rsid (character),
#     a0 (other/ref allele), a1 (effect allele), freq (allele-1 frequency).
#   pops: character vector of super-population codes to return (default all 6).
# Returns:
#   named numeric vector of weights summing to 1, one entry per pop in `pops`.
#   NA if AF-projection cannot be computed (e.g. too few matched SNPs).

suppressPackageStartupMessages({
  library(data.table)
  library(bigsnpr)
  library(bigreadr)
})

BIGSNPR_DIR <- '/users/k1806347/oliverpainfel/Data/bigsnpr'

# Privé's ancestry-vignette PC-correction constant (from bigsnpr docs).
.AFPROJ_CORRECTION <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                        1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

# Fine → coarse regional groups (folds Scandinavia/UK/Ireland → "Europe (North West)",
# Europe South East/North East → "Europe (East)"). Matches Round 7 usage.
.coarse_group <- function(fine) {
  g <- fine
  g[g %in% c('Scandinavia','United Kingdom','Ireland')]     <- 'Europe (North West)'
  g[g %in% c('Europe (South East)','Europe (North East)')]  <- 'Europe (East)'
  g
}

# Coarse group → GenoPred super-population code.
.SUPER_MAP <- c(
  'Africa (West)' = 'AFR', 'Africa (South)' = 'AFR',
  'Africa (East)' = 'AFR', 'Africa (North)' = 'AFR',
  'Middle East' = 'MID', 'Ashkenazi' = 'EUR', 'Italy' = 'EUR',
  'Finland' = 'EUR',
  'Europe (East)' = 'EUR', 'Europe (North West)' = 'EUR',
  'Europe (South West)' = 'EUR',
  'South America' = 'AMR',
  'Sri Lanka' = 'CSA', 'Pakistan' = 'CSA', 'Bangladesh' = 'CSA',
  'Asia (East)' = 'EAS', 'Japan' = 'EAS', 'Philippines' = 'EAS'
)

compute_afproj_weights <- function(gwas_freq,
                                   pops = c('EUR','EAS','AFR','CSA','AMR','MID'),
                                   match_min_prop = 0.05,
                                   redistribute_mid_to = NA_character_) {
  stopifnot(is.data.table(gwas_freq),
            all(c('chr','pos','rsid','a0','a1','freq') %in% names(gwas_freq)))

  all_freq   <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'ref_freqs.csv.gz'))
  projection <- bigreadr::fread2(file.path(BIGSNPR_DIR, 'projection.csv.gz'))

  # snp_match wants a beta column (any dummy).
  gf <- copy(gwas_freq)
  gf[, beta := 1]
  # Match on rsid rather than chr+pos so this works across genome builds
  # (GBMI is b38 while all_freq is b37).
  matched <- bigsnpr::snp_match(as.data.frame(gf), all_freq[, 1:5],
                                match.min.prop = match_min_prop,
                                join_by_pos    = FALSE)
  if (nrow(matched) < 100) {
    warning(sprintf('AF-projection: only %d SNPs matched; returning NA', nrow(matched)))
    out <- setNames(rep(NA_real_, length(pops)), pops)
    return(out)
  }
  matched$freq <- ifelse(matched$beta < 0, 1 - matched$freq, matched$freq)

  res <- bigsnpr::snp_ancestry_summary(
    freq          = matched$freq,
    info_freq_ref = all_freq[matched$`_NUM_ID_`, -(1:5)],
    projection    = projection[matched$`_NUM_ID_`, -(1:5)],
    correction    = .AFPROJ_CORRECTION
  )

  fine_pops <- colnames(all_freq)[-(1:5)]
  grp_fct   <- factor(.coarse_group(fine_pops),
                      levels = unique(.coarse_group(fine_pops)))
  by_coarse <- tapply(res, grp_fct, sum)
  by_super  <- tapply(by_coarse, .SUPER_MAP[names(by_coarse)], sum, default = 0)

  w <- setNames(numeric(length(pops)), pops)
  for (s in pops) if (s %in% names(by_super)) w[[s]] <- as.numeric(by_super[[s]])

  # Optional: redistribute MID mass into a target super-pop if MID is not a
  # candidate. Set redistribute_mid_to='EUR' to mirror the R7 behaviour.
  if (!is.na(redistribute_mid_to) && 'MID' %in% names(w) && w[['MID']] > 0) {
    if (!(redistribute_mid_to %in% names(w)))
      stop(sprintf("redistribute_mid_to='%s' not in `pops`", redistribute_mid_to))
    w[[redistribute_mid_to]] <- w[[redistribute_mid_to]] + w[['MID']]
    w[['MID']] <- 0
  }

  s <- sum(w); if (s > 0) w / s else w
}

# convenience: two composition scalars per prompt §4
composition_scalars <- function(weights, per_pop_af = NULL) {
  # weights: named numeric summing to 1 (may contain zeros).
  # per_pop_af: named list of numeric vectors of AF per pop over shared SNPs;
  #             all vectors must have the same length. If NULL, only the effective
  #             ancestry-groups scalar is returned.
  w <- weights[weights > 0]
  eff_groups <- 1 / sum(w^2)

  mean_B_diag <- NA_real_
  if (!is.null(per_pop_af)) {
    common <- intersect(names(per_pop_af), names(w))
    if (length(common) >= 2) {
      F <- do.call(cbind, per_pop_af[common])
      wv <- w[common] / sum(w[common])
      fbar <- as.vector(F %*% wv)
      dF <- F - fbar
      # B(i,i) = 4 Σ_p w_p (f_p,i − f̄_i)^2  — mean over SNPs of the diagonal
      diag_B <- 4 * as.vector((dF^2) %*% wv)
      mean_B_diag <- mean(diag_B, na.rm = TRUE)
    }
  }
  c(eff_groups = eff_groups, mean_B_diag = mean_B_diag)
}
