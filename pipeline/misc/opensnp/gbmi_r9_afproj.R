#!/usr/bin/env Rscript
# gbmi_r9_afproj.R — Round 9 Section 3.
#
# Re-run the ancestry-composition estimation against UKB per-population
# reference allele frequencies (as recorded in the UKB block RDS f_{POP}
# fields), and report whether the AF-projection weights shift from Round 8.
#
# Round 8 used bigsnpr::snp_ancestry_summary against Privé's 1KG+HGDP
# fine-population reference (16 populations, projection loadings). UKB
# provides only 6 super-populations without projection loadings, so a
# direct swap into snp_ancestry_summary is not possible. Instead:
# fit non-negative weights via a constrained quadratic program
#   min ||f_meta - F * w||^2   s.t.  w >= 0, sum(w) = 1
# where F is the (n_SNPs × 6) matrix of per-population reference AFs.
#
# Two runs per trait:
#   w_qp_r8  — QP against Round 8's 1KG+HGDP per-pop AFs (block RDS f_{POP})
#   w_qp_r9  — QP against UKB per-pop AFs (block RDS f_{POP})
# and comparison against gbmi_r8_weights.csv's w_afproj (the snp_ancestry_summary
# result — different estimator; only for context).
#
# Emits gbmi_r9_weights.csv with columns
#   trait, pop, w_qp_r8, w_qp_r9, w_r8_afproj_stored, delta_qp = w_qp_r9 - w_qp_r8

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(quadprog)
})

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
source(file.path(MISC, 'meld_paths.R'))
OUT_DIR <- r_results('r9')
R8_LD_DIR <- file.path(meld_ld_1kg_hgdp(), 'chr22')
R9_LD_DIR <- file.path(meld_ld_ukb(), 'chr22')

# Pops in the QP:
POPS_R8 <- c('EUR','EAS','AFR','CSA','AMR')      # Round 8 blocks have no MID
POPS_R9 <- c('EUR','EAS','AFR','CSA','AMR','MID') # UKB blocks have MID

# ---------------------------------------------------------------------------
load_freq_from_blocks <- function(block_dir, pops) {
  files <- sort(list.files(block_dir, pattern = '^block_.*\\.rds$', full.names = TRUE))
  frames <- lapply(files, function(fp) {
    bl <- readRDS(fp)
    dt <- data.table(rsid = bl$SNP)
    for (P in pops) dt[[P]] <- bl[[sprintf('f_%s', P)]]
    dt
  })
  d <- rbindlist(frames)
  # For duplicated rsids across blocks, keep the first (they should agree)
  d[!duplicated(d$rsid)]
}

qp_fit <- function(f_query, F_ref) {
  # min 0.5 * w' * (2 F'F) w - 2 (F' f_query)' w, s.t. sum(w)=1, w >= 0
  # solve.QP form: min 0.5 w'Dw - d'w  with A'w >= b, first meq are equalities.
  D <- 2 * crossprod(F_ref)  # p x p
  d <- 2 * crossprod(F_ref, f_query)  # p
  p <- ncol(F_ref)
  # constraints: sum(w) = 1 (equality), w >= 0 (inequalities)
  A <- cbind(rep(1, p), diag(p))
  b <- c(1, rep(0, p))
  # small ridge for numerical stability
  D <- D + diag(1e-6, p, p)
  sol <- solve.QP(D, d, A, b, meq = 1)
  w <- sol$solution
  w[w < 0] <- 0
  w / sum(w)
}

# ---------------------------------------------------------------------------
cat('=== loading reference AFs ===\n')
r8_af <- load_freq_from_blocks(R8_LD_DIR, POPS_R8)
r9_af <- load_freq_from_blocks(R9_LD_DIR, POPS_R9)
cat(sprintf('  R8: %d unique rsids across 5 pops\n', nrow(r8_af)))
cat(sprintf('  R9: %d unique rsids across 6 pops\n', nrow(r9_af)))

# ---------------------------------------------------------------------------
TRAITS <- c('Asthma','COPD','Gout','HF','IPF','Stroke','VTE')
r8_stored <- fread(file.path(r_results('r8'), 'gbmi_r8_weights.csv'))

rows <- list()
for (TR in TRAITS) {
  d <- readRDS(r8_harmonised(TR))
  cat(sprintf('\n=== %s (n=%d) ===\n', TR, nrow(d)))

  # QP against R8 reference
  d_r8 <- merge(d[, .(rsid, meta_af)], r8_af, by = 'rsid')
  F8 <- as.matrix(d_r8[, ..POPS_R8])
  w8 <- qp_fit(d_r8$meta_af, F8)
  names(w8) <- POPS_R8
  cat('QP vs R8 1KG+HGDP AFs:\n'); print(round(w8, 4))

  # QP against R9 (UKB) reference
  d_r9 <- merge(d[, .(rsid, meta_af)], r9_af, by = 'rsid')
  F9 <- as.matrix(d_r9[, ..POPS_R9])
  w9 <- qp_fit(d_r9$meta_af, F9)
  names(w9) <- POPS_R9
  cat('QP vs R9 UKB AFs:\n'); print(round(w9, 4))

  # Stored R8 afproj from gbmi_r8_weights.csv (snp_ancestry_summary, different estimator)
  stored <- r8_stored[trait == TR]
  w_stored <- setNames(rep(0, length(POPS_R9)), POPS_R9)
  for (P in POPS_R9) {
    v <- stored[pop == P, w_afproj]
    if (length(v) == 1L && is.finite(v)) w_stored[P] <- v
  }

  # Assemble a per-pop row
  for (P in POPS_R9) {
    rows[[length(rows) + 1L]] <- data.table(
      trait = TR, pop = P,
      w_qp_r8 = if (P %in% POPS_R8) w8[[P]] else NA_real_,
      w_qp_r9 = if (P %in% POPS_R9) w9[[P]] else NA_real_,
      w_r8_afproj_stored = w_stored[[P]]
    )
  }
}

W <- rbindlist(rows)
W[, delta_qp := w_qp_r9 - w_qp_r8]

# Report movement
cat('\n\n=== movement in QP weights when swapping reference AFs (R9 UKB - R8 1KG+HGDP) ===\n')
mv <- dcast(W, trait ~ pop, value.var = 'delta_qp')
print(mv)

cat('\n=== max |delta| per trait ===\n')
mx <- W[!is.na(delta_qp), .(max_abs_delta = max(abs(delta_qp), na.rm = TRUE)), by = trait]
print(mx)

# The MID column in delta_qp is NA (no R8 QP for MID). Report R9 MID mass separately.
cat('\n=== MID mass under R9 (was folded into EUR in R8) ===\n')
print(W[pop == 'MID', .(trait, w_qp_r9)])

out <- file.path(OUT_DIR, 'gbmi_r9_weights.csv')
fwrite(W, out)
cat(sprintf('\nwrote %s\n', out))
