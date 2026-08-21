#!/usr/bin/env Rscript
# Round 3: construct synthetic mixed-ancestry GWAS from Yengo EUR + AFR
# per-population sumstats, at a set of ancestry proportions p.
#
# Per SNP, after allele harmonisation between the two source files:
#   AF_mix = p_eur * AF_eur + p_afr * AF_afr
#   b_mix  = p_eur * b_eur  + p_afr * b_afr
#   N_mix  = fixed (549000 = median Yengo-EUR N) so all mixtures share
#           the same "cohort size", isolating the ancestry effect from
#           any N-scaling of power
#   SE_mix = 1 / sqrt(N_mix * 2 * AF_mix * (1 - AF_mix))     [sample-size formula]
#   P_mix  = 2 * pnorm(-abs(b_mix / SE_mix))

suppressPackageStartupMessages({
  library(data.table)
})

sumstats_dir <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/opensnp_test'
eur <- fread(file.path(sumstats_dir, 'yengo_2022_height_eur.txt'))
afr <- fread(file.path(sumstats_dir, 'yengo_2022_height_afr.txt'))

# Drop rows with missing rsid and de-duplicate — the Yengo files have a small
# number of NA variant_id rows plus a few genuine duplicates (~560 per file).
# A merge on variant_id would produce spurious many-to-many joins on the NAs.
eur <- eur[!is.na(variant_id)][!duplicated(variant_id)]
afr <- afr[!is.na(variant_id)][!duplicated(variant_id)]
cat(sprintf('After NA-drop + dedup: EUR %d, AFR %d\n', nrow(eur), nrow(afr)))
h <- merge(
  eur[, .(variant_id, chromosome, base_pair_location,
          ea_eur = effect_allele, oa_eur = other_allele,
          beta_eur = beta, freq_eur = effect_allele_frequency,
          se_eur = standard_error, p_eur_in = p_value, n_eur = n)],
  afr[, .(variant_id,
          ea_afr = effect_allele, oa_afr = other_allele,
          beta_afr = beta, freq_afr = effect_allele_frequency,
          se_afr = standard_error, p_afr_in = p_value, n_afr = n)],
  by = 'variant_id'
)

cat(sprintf('Merged %d SNPs (EUR %d, AFR %d)\n',
            nrow(h), nrow(eur), nrow(afr)))

# Same-allele SNPs → keep as-is
same <- h$ea_eur == h$ea_afr & h$oa_eur == h$oa_afr
# Flipped-allele SNPs (ea_eur == oa_afr and oa_eur == ea_afr) → flip AFR
flip <- h$ea_eur == h$oa_afr & h$oa_eur == h$ea_afr
# Anything else (ambiguous / different strand / mismatched) → drop
drop <- !(same | flip)

cat(sprintf('  same-alleles: %d\n  flipped:      %d\n  dropped:      %d\n',
            sum(same), sum(flip), sum(drop)))

h <- h[!drop]
h[flip[!drop], `:=`(
  beta_afr = -beta_afr,
  freq_afr = 1 - freq_afr
)]

# Reference N for all mixtures (isolates ancestry effect from power effect)
N_MIX <- 549000L

# Build a mixture at proportion p_eur (p_afr = 1 - p_eur)
build_mix <- function(p_eur) {
  p_afr <- 1 - p_eur
  af    <- p_eur * h$freq_eur + p_afr * h$freq_afr
  b     <- p_eur * h$beta_eur + p_afr * h$beta_afr
  se    <- 1 / sqrt(N_MIX * 2 * af * (1 - af))
  P     <- 2 * pnorm(-abs(b / se))
  data.table(
    chromosome              = h$chromosome,
    base_pair_location      = h$base_pair_location,
    effect_allele           = h$ea_eur,
    other_allele            = h$oa_eur,
    beta                    = b,
    standard_error          = se,
    effect_allele_frequency = af,
    p_value                 = P,
    variant_id              = h$variant_id,
    n                       = N_MIX
  )
}

mix_points <- list(
  eur100 = 1.00,
  eur75  = 0.75,
  eur50  = 0.50,
  eur25  = 0.25,
  eur00  = 0.00
)

for (nm in names(mix_points)) {
  p_eur <- mix_points[[nm]]
  out <- build_mix(p_eur)
  # Drop SNPs with degenerate AF (0 or 1) — no information, would divide by zero in SE
  keep <- out$effect_allele_frequency > 0 & out$effect_allele_frequency < 1
  out  <- out[keep]

  path <- file.path(sumstats_dir, sprintf('yengo_2022_height_mix_%s.txt', nm))
  fwrite(out, path, sep = '\t')
  cat(sprintf('wrote %s (%d SNPs, p_eur=%.2f)\n', path, nrow(out), p_eur))
}
