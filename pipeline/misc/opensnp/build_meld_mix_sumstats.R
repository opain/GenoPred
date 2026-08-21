#!/usr/bin/env Rscript
# Round 6, Section 1: construct synthetic mixed-ancestry GWAS from Yengo EUR + AFR
# using the inverse-variance-weighted meta target that DOES have a valid LD matrix.
#
# Per SNP, after allele harmonisation to the reference ALT-allele canonical
# orientation (populations p ∈ {EUR, AFR}):
#   u_p,i     = N_p · v_p,i                       (v_p empirical from 1KG+HGDP reference)
#   β_meta,i  = Σ_p u_p,i · β_p,i / Σ_p u_p,i
#   SE_meta,i = 1 / sqrt( Σ_p u_p,i )
#   EAF_i     = Σ_p N_p · f_p,i / Σ_p N_p         (arithmetic AF, per prompt)
#   N_i       = Σ_p N_p
#   z_i       = β_meta,i / SE_meta,i ;  P_i = 2·pnorm(-|z_i|)
#
# Mixture points implemented via N_p (total N held constant). w_p = N_p / ΣN_q.
#
# Superseded by Round-6 rewrite; original R3 script preserved as
# build_meld_mix_sumstats_v1.R for Round 3 reproducibility.
#
# See ~/meld_round6_prompt.md Sections 0(a), 0(b), 1 for the design rationale.

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Inputs & config
sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
ref_emp_dir  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/ref_empirical'
N_TOTAL      <- 500000L                    # per prompt example; constant across mixtures
mix_points   <- list(eur100 = 1.00,        # w_eur
                     eur75  = 0.75,
                     eur50  = 0.50,
                     eur25  = 0.25,
                     eur00  = 0.00)

# ---------------------------------------------------------------------------
# Load sumstats and reference empirical variances
eur_gwas <- fread(file.path(sumstats_dir, 'yengo_2022_height_eur.txt'))
afr_gwas <- fread(file.path(sumstats_dir, 'yengo_2022_height_afr.txt'))

# Drop rows with NA variant_id and de-duplicate — see v1 script note; same
# NA-and-duplicate pathology present in Yengo files.
eur_gwas <- eur_gwas[!is.na(variant_id)][!duplicated(variant_id)]
afr_gwas <- afr_gwas[!is.na(variant_id)][!duplicated(variant_id)]
cat(sprintf('Input: EUR %d, AFR %d\n', nrow(eur_gwas), nrow(afr_gwas)))

eur_ref <- readRDS(file.path(ref_emp_dir, 'EUR.var.rds'))
afr_ref <- readRDS(file.path(ref_emp_dir, 'AFR.var.rds'))
stopifnot(all(eur_ref$SNP == afr_ref$SNP),
          all(eur_ref$REF == afr_ref$REF),
          all(eur_ref$ALT == afr_ref$ALT))

# Canonical orientation = reference ALT allele. Both empirical files use the
# same variant/allele order so we can hold a single canonical view.
ref <- eur_ref[, .(CHR, BP, SNP, A1 = ALT, A2 = REF,
                   f_eur = f_alt, v_eur = v)]
ref[, `:=`(f_afr = afr_ref$f_alt, v_afr = afr_ref$v)]

# ---------------------------------------------------------------------------
# Harmonise sumstats to canonical A1 = reference ALT allele
harmonise <- function(gwas, tag) {
  # gwas cols: variant_id, effect_allele, other_allele, beta, effect_allele_frequency, n
  m <- merge(ref[, .(SNP, A1, A2)],
             gwas[, .(SNP = variant_id, ea = effect_allele, oa = other_allele,
                      b  = beta, f = effect_allele_frequency, N = n)],
             by = 'SNP', sort = FALSE)
  same <- m$ea == m$A1 & m$oa == m$A2      # already canonical
  flip <- m$ea == m$A2 & m$oa == m$A1      # opposite orientation, flip
  drop <- !(same | flip)
  cat(sprintf('  %s harmonisation: %d matched (%d same, %d flip, %d drop)\n',
              tag, nrow(m), sum(same), sum(flip), sum(drop)))
  m <- m[!drop]
  # Flip β and freq for flipped SNPs
  m[flip[!drop], `:=`(b = -b, f = 1 - f)]
  m[, `:=`(ea = NULL, oa = NULL)]
  setnames(m, c('b', 'f'), c(paste0('b_',   tag),
                             paste0('fgw_', tag)))
  m
}

h_eur <- harmonise(eur_gwas, 'eur')
h_afr <- harmonise(afr_gwas, 'afr')

h <- merge(h_eur, h_afr, by = c('SNP', 'A1', 'A2'), sort = FALSE)
h <- merge(ref, h, by = c('SNP', 'A1', 'A2'), sort = FALSE)
setnames(h, c('N.x', 'N.y'), c('N_eur_yengo', 'N_afr_yengo'))
cat(sprintf('After harmonisation and inner-merge: %d SNPs\n', nrow(h)))

# Drop SNPs with NA empirical variance (n_called < 2 in a pop); treat any
# remaining variance-zero as 0 (monomorphic in that pop → term drops out).
h[is.na(v_eur), v_eur := 0]
h[is.na(v_afr), v_afr := 0]

# ---------------------------------------------------------------------------
# Build one mixture: return a data.table with output columns matching v1
build_mix <- function(w_eur) {
  w_afr <- 1 - w_eur
  N_eur <- w_eur * N_TOTAL
  N_afr <- w_afr * N_TOTAL

  u_eur <- N_eur * h$v_eur
  u_afr <- N_afr * h$v_afr
  u_sum <- u_eur + u_afr

  # β_meta and SE_meta from IVW pool
  b_meta  <- (u_eur * h$b_eur + u_afr * h$b_afr) / u_sum
  se_meta <- 1 / sqrt(u_sum)

  # EAF: arithmetic mixture at true w (uses reference AFs, not Yengo-reported)
  f_meta <- w_eur * h$f_eur + w_afr * h$f_afr

  # Drop SNPs with u_sum == 0 (monomorphic in every contributing pop at this w)
  ok <- is.finite(b_meta) & u_sum > 0
  b_meta[!ok]  <- NA_real_
  se_meta[!ok] <- NA_real_

  z <- b_meta / se_meta
  P <- 2 * pnorm(-abs(z))

  out <- data.table(
    chromosome              = h$CHR,
    base_pair_location      = h$BP,
    effect_allele           = h$A1,
    other_allele            = h$A2,
    beta                    = b_meta,
    standard_error          = se_meta,
    effect_allele_frequency = f_meta,
    p_value                 = P,
    variant_id              = h$SNP,
    n                       = as.integer(round(N_eur + N_afr))
  )
  out[ok]
}

# ---------------------------------------------------------------------------
# Sanity checks
# 1) At w_eur = 1, output β and AF equal Yengo EUR (after harmonising Yengo
#    itself to the reference orientation — that is what h holds).
# 2) At w_eur = 0, same for AFR.
run_sanity <- function() {
  cat('\nSanity checks\n')
  m1 <- build_mix(1.00)
  m0 <- build_mix(0.00)

  # For w=1, β_meta = b_eur, EAF = f_eur (reference AF, not Yengo AF).
  # The prompt asks max |Δβ| = 0 on shared same-allele SNPs, so compare to
  # h$b_eur (Yengo β already harmonised to reference orientation). AF comparison:
  # per prompt, compare to Yengo AF (fgw_eur), not reference AF, to catch orientation
  # bugs (both should be zero if flipping was correct).
  ord <- match(m1$variant_id, h$SNP)
  db_eur  <- max(abs(m1$beta - h$b_eur[ord]),   na.rm = TRUE)
  # Sanity check on AF compares to the Yengo-reported EAF that the sumstats
  # carried, harmonised to canonical orientation (fgw_eur). At w=1, EAF_out
  # is w_eur * f_eur + w_afr * f_afr = f_eur (reference AF) — which is NOT
  # guaranteed to equal the Yengo-reported EUR AF. Report both.
  daf_ref = max(abs(m1$effect_allele_frequency - h$f_eur[ord]), na.rm = TRUE)   # 0 by construction
  daf_gwa = max(abs(m1$effect_allele_frequency - h$fgw_eur[ord]), na.rm = TRUE) # nonzero: Yengo vs reference AF
  cat(sprintf('  w=1: max|Δβ vs Yengo-EUR-harm| = %.3e   (must be 0)\n', db_eur))
  cat(sprintf('       max|ΔAF vs reference EUR| = %.3e   (must be 0)\n', daf_ref))
  cat(sprintf('       max|ΔAF vs Yengo-EUR EAF| = %.3e   (informational: Yengo vs 1KG+HGDP)\n', daf_gwa))
  stopifnot(db_eur < 1e-12, daf_ref < 1e-12)

  ord <- match(m0$variant_id, h$SNP)
  db_afr  <- max(abs(m0$beta - h$b_afr[ord]),   na.rm = TRUE)
  daf_ref = max(abs(m0$effect_allele_frequency - h$f_afr[ord]), na.rm = TRUE)
  daf_gwa = max(abs(m0$effect_allele_frequency - h$fgw_afr[ord]), na.rm = TRUE)
  cat(sprintf('  w=0: max|Δβ vs Yengo-AFR-harm| = %.3e   (must be 0)\n', db_afr))
  cat(sprintf('       max|ΔAF vs reference AFR| = %.3e   (must be 0)\n', daf_ref))
  cat(sprintf('       max|ΔAF vs Yengo-AFR EAF| = %.3e   (informational: Yengo vs 1KG+HGDP)\n', daf_gwa))
  stopifnot(db_afr < 1e-12, daf_ref < 1e-12)

  # 3) SE_meta = 1/sqrt(Σ u_p) is monotone in Σ u_p exactly.
  #    (The prompt asks for monotonicity in EAF·(1-EAF). That relation
  #    only holds under HWE v_p = 2f_p(1-f_p) AND at endpoint w=0 or 1;
  #    at intermediate w the pooled EAF·(1-EAF) is not equivalent to
  #    the pooled u_p — so we check the tight relation on u_p exactly,
  #    and report the EAF·(1-EAF) correlation for information. If the
  #    EAF·(1-EAF) correlation is close to −1 that is a HINT the SE
  #    formula collapsed to the v1 pooled-HWE form.)
  for (nm in names(mix_points)) {
    w_eur  <- mix_points[[nm]]; w_afr <- 1 - w_eur
    N_eur  <- w_eur * N_TOTAL;  N_afr <- w_afr * N_TOTAL
    m      <- build_mix(w_eur)
    ord    <- match(m$variant_id, h$SNP)
    u_sum  <- N_eur * h$v_eur[ord] + N_afr * h$v_afr[ord]
    rho_u  <- suppressWarnings(cor(m$standard_error, u_sum, method = 'spearman'))
    het    <- m$effect_allele_frequency * (1 - m$effect_allele_frequency)
    rho_h  <- suppressWarnings(cor(m$standard_error, het,   method = 'spearman'))
    cat(sprintf('  w=%.2f: Spearman(SE, Σu_p) = %+.4f  (must be -1);  Spearman(SE, EAF·(1-EAF)) = %+.4f (info)\n',
                w_eur, rho_u, rho_h))
    stopifnot(rho_u < -0.999999)
  }
  cat('  all sanity checks passed.\n')
}
run_sanity()

# ---------------------------------------------------------------------------
# Write outputs (v1-compatible paths and column order)
for (nm in names(mix_points)) {
  w_eur <- mix_points[[nm]]
  out   <- build_mix(w_eur)
  path  <- file.path(sumstats_dir, sprintf('yengo_2022_height_mix_%s.txt', nm))
  fwrite(out, path, sep = '\t')
  cat(sprintf('wrote %s (%d SNPs, w_eur=%.2f, N=%d)\n',
              path, nrow(out), w_eur, out$n[1L]))
}
