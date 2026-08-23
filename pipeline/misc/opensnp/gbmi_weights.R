#!/usr/bin/env Rscript
# gbmi_weights.R — Round 8 Section 4.
#
# Per trait: compute mixture weights two ways:
#   (a) AF-projection: bigsnpr::snp_ancestry_summary on meta_af against the
#       1KG+HGDP fine-pop reference, rolled up to super-pops.
#   (b) Reported per-ancestry N: median per-SNP N per arm, normalised to 1.
# Plus two composition scalars:
#   - Effective number of ancestry groups: 1 / Σ w_p²
#   - Mean between-population variance: mean over SNPs of 4·Σ_p w_p (f_p,i − f̄_i)²
# For traits whose meta includes a population with no arm (SAS in 6/7 traits,
# MID in Asthma), per-arm AF is unavailable — fall back on reference-panel AF.
#
# Writes:
#   gbmi_r8_weights.csv        one row per (trait, pop, source)
#   gbmi_r8_composition.csv    one row per trait (eff_groups, mean_B_diag,
#                              plus the trait's weight vector for reference)

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages(library(data.table))

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r8')
REF_FREQ_DIR <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref/freq_files'
source(file.path(MISC, 'af_projection.R'))

TRAITS <- c('Asthma','COPD','Gout','HF','IPF','Stroke','VTE')
# For each trait, which super-pops are in the meta (from manifest §4).
META_POPS <- list(
  Asthma = c('EUR','EAS','AFR','CSA','AMR','MID'),
  COPD   = c('EUR','EAS','AFR','CSA','AMR'),
  Gout   = c('EUR','EAS','AFR','CSA','AMR'),
  HF     = c('EUR','EAS','AFR','CSA','AMR'),
  IPF    = c('EUR','EAS','AFR','CSA','AMR'),
  Stroke = c('EUR','EAS','AFR','CSA','AMR'),
  VTE    = c('EUR','EAS','AFR','CSA','AMR')
)
# Per-trait: which pops have arms (used for reported-N weights and per-pop AF).
ARM_POPS <- list(
  Asthma = c('EUR','EAS','AFR','CSA','AMR'),   # SAS→CSA; no MID arm
  COPD   = c('EUR','EAS','AFR','AMR'),          # no SAS arm
  Gout   = c('EUR','EAS','AFR','AMR'),
  HF     = c('EUR','EAS','AFR','AMR'),
  IPF    = c('EUR','EAS','AFR','AMR'),
  Stroke = c('EUR','EAS','AFR','AMR'),
  VTE    = c('EUR','EAS','AFR','AMR')
)

# Reference-panel AF fallback per pop (chr22 only, on the same rsid set).
load_ref_af_chr22 <- function(pop, rsids) {
  fp <- file.path(REF_FREQ_DIR, pop, sprintf('ref.%s.chr22.afreq', pop))
  af <- fread(fp, showProgress = FALSE)
  # columns: #CHROM ID REF ALT ALT_FREQS OBS_CT
  setnames(af, '#CHROM', 'chr')
  af[match(rsids, af$ID), .(rsid = rsids, ref_af = REF, alt_af = ALT, af = ALT_FREQS)]
}

# ---------------------------------------------------------------------------
weights_all      <- list()
composition_all  <- list()
per_pop_af_used  <- list()   # for cross-checks

for (TRAIT in TRAITS) {
  cat(sprintf('\n=== %s ===\n', TRAIT))
  rds <- r8_harmonised(TRAIT)
  d <- readRDS(rds)
  cat(sprintf('n variants: %d\n', nrow(d)))

  meta_pops <- META_POPS[[TRAIT]]
  arm_pops  <- ARM_POPS [[TRAIT]]

  # --- Route (a): AF-projection weights -----------------------------------
  gwas_freq <- d[, .(chr, pos, rsid, a0 = ref, a1 = alt, freq = meta_af)]
  w_afproj <- compute_afproj_weights(gwas_freq, pops = c('EUR','EAS','AFR','CSA','AMR','MID'))
  cat('AF-projection weights (super-pop):\n')
  print(round(w_afproj, 4))

  # --- Route (b): Reported per-ancestry N ---------------------------------
  N_by_arm <- setNames(numeric(length(arm_pops)), arm_pops)
  for (P in arm_pops) {
    col <- sprintf('N_%s', P)
    if (!(col %in% names(d))) stop(sprintf('missing column %s', col))
    N_by_arm[P] <- as.numeric(median(d[[col]], na.rm = TRUE))
  }
  w_N <- N_by_arm / sum(N_by_arm)
  cat('Reported-N weights (arm pops):\n')
  print(round(w_N, 4))

  # --- Save side-by-side into weights_all ---------------------------------
  # For AF-projection, all 6 super-pops potentially non-zero. For reported-N,
  # only arm_pops. Fill NA for the pops that don't apply to the other route.
  all_pops <- union(names(w_afproj), names(w_N))
  weights_all[[TRAIT]] <- data.table(
    trait = TRAIT,
    pop = all_pops,
    w_afproj = unname(w_afproj[all_pops]),
    w_N      = unname(w_N     [all_pops]),
    N_median = unname(N_by_arm[all_pops])
  )
  # w_N is NA for pops without an arm; N_median is NA there too.
  weights_all[[TRAIT]][, discrepancy_afproj_minus_N :=
                        ifelse(is.na(w_N), NA_real_, w_afproj - w_N)]

  # --- Composition scalar (i): effective ancestry groups -----------------
  eff_groups <- 1 / sum(w_afproj[w_afproj > 0]^2)

  # --- Composition scalar (ii): mean between-pop variance ----------------
  # Build a per-pop AF matrix over the trait's meta pops. For pops with an arm,
  # use the harmonised arm AF. For pops without an arm (SAS in non-Asthma, MID
  # in Asthma), use reference-panel AF matched by rsid.
  per_pop_af <- list()
  for (P in meta_pops) {
    if (P %in% arm_pops) {
      col <- sprintf('af_%s', P)
      per_pop_af[[P]] <- d[[col]]
    } else {
      # Reference-panel AF for this pop, aligned to d's rsid + REF/ALT.
      ra <- load_ref_af_chr22(P, d$rsid)
      # Refpop AF may or may not match d's REF/ALT orientation.
      # d has been aligned to HapMap3 ref/alt; the ref freq file is on the same
      # HM3 panel (built with the same pvar), so orientation is identical.
      # Cross-check: same REF and ALT.
      if (any(ra$ref_af != d$ref | ra$alt_af != d$alt, na.rm = TRUE))
        stop(sprintf('ref-panel AF orientation differs from harmonised orientation for %s', P))
      per_pop_af[[P]] <- ra$af
      cat(sprintf('  used reference-panel AF for %s (%d SNPs; NA=%d)\n',
                  P, length(ra$af), sum(is.na(ra$af))))
    }
  }
  # Align the pop weights to the meta_pops set. w_afproj covers all 6, restrict.
  w_here <- w_afproj[meta_pops]
  w_here <- w_here / sum(w_here)   # renormalise inside the trait's meta

  cs <- composition_scalars(w_here, per_pop_af)
  composition_all[[TRAIT]] <- data.table(
    trait          = TRAIT,
    n_variants     = nrow(d),
    eff_groups     = as.numeric(cs[['eff_groups']]),
    mean_B_diag    = as.numeric(cs[['mean_B_diag']]),
    top_pop        = names(which.max(w_here)),
    top_pop_weight = max(w_here)
  )

  # Remember per-pop AF used (for downstream cross-checks and the notebook).
  per_pop_af_used[[TRAIT]] <- list(pops = meta_pops, af = per_pop_af,
                                   w = as.numeric(w_here))
  cat(sprintf('eff_groups = %.3f, mean_B_diag = %.6f\n',
              cs[['eff_groups']], cs[['mean_B_diag']]))
}

W  <- rbindlist(weights_all)
CS <- rbindlist(composition_all)

fwrite(W,  file.path(OUT_DIR, 'gbmi_r8_weights.csv'))
fwrite(CS, file.path(OUT_DIR, 'gbmi_r8_composition.csv'))
saveRDS(per_pop_af_used, file.path(MELD_DATA, 'gbmi_r8_harmonised', 'gbmi_r8_per_pop_af.rds'))

cat('\n\n=== Weights summary (agreement AF-proj vs reported-N) ===\n')
# Focus on arm-covered pops per trait (where both routes apply).
W_arm <- W[!is.na(w_N)]
W_arm[, disc_txt := sprintf('%+0.4f', w_afproj - w_N)]
print(dcast(W_arm, trait ~ pop, value.var = 'disc_txt', fill = ''))
cat('\nMax |discrepancy| per trait:\n')
print(W_arm[, .(max_abs_disc = max(abs(w_afproj - w_N), na.rm = TRUE)), by = trait])

cat('\n=== Composition scalars ===\n')
print(CS)
cat('\nDONE\n')
