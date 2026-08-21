#!/usr/bin/env Rscript
# Round 7 Section 1 helper — intersect variant sets across candidate PGS panels
# so any downstream comparison (M5 weight-vector agreement, or target R²) runs
# on an identical SNP set. Reports intersection size and per-panel drop counts.
#
# Inputs are per-panel PGS score files as produced by GenoPred (e.g.
# `.raw.profiles` from ptclump/dbslmm/quickprs; the "score" file with per-SNP
# weights per model). The script intersects on SNP ID.
#
# Usage:
#   Rscript intersect_pgs_variants.R <config.tsv>
# where <config.tsv> is a tab-separated 3-column file:
#   panel_label   score_file_path   effect_col
# effect_col is the column name holding the per-SNP effect (BETA / SCORE etc.).
#
# Writes:
#   <outdir>/intersect_snplist.txt         one SNP per line (intersection)
#   <outdir>/intersect_summary.csv         per-panel drop counts, intersection size
#   <outdir>/intersect_weights.rds         list of per-panel weights on the intersection
#
# The idea is that any subsequent weight-vector correlation, target-R² etc. can
# be computed against `intersect_weights.rds` directly and is guaranteed to run
# on the same SNP set for every panel.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop('usage: intersect_pgs_variants.R <config.tsv> [<outdir>]')
cfg_path <- args[[1L]]
out_dir  <- if (length(args) >= 2L) args[[2L]] else dirname(cfg_path)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cfg <- fread(cfg_path, header = TRUE)
required <- c('panel_label', 'score_file_path', 'effect_col')
if (!all(required %in% names(cfg))) {
  stop('config must have columns: ', paste(required, collapse = ', '))
}

# Read each panel's per-SNP weights.
per_panel <- list()
for (i in seq_len(nrow(cfg))) {
  lbl <- cfg$panel_label[[i]]
  fp  <- cfg$score_file_path[[i]]
  ec  <- cfg$effect_col[[i]]
  if (!file.exists(fp)) stop('score file not found: ', fp)
  dt <- fread(fp)
  snp_col <- intersect(c('SNP', 'ID', 'variant_id', 'MarkerName', 'RSID'), names(dt))
  if (length(snp_col) == 0L) stop('cannot find SNP identifier column in ', fp,
                                  ' (looked for: SNP, ID, variant_id, MarkerName, RSID)')
  snp_col <- snp_col[[1L]]
  if (!(ec %in% names(dt))) stop('effect column ', ec, ' missing from ', fp)
  # Deduplicate SNP rows — some methods produce multiple rows per SNP (thresholds/h²
  # variants). Aggregate by SNP using the first available value.
  dt <- dt[, .(effect = get(ec)[1L]), by = c('SNP' = snp_col)]
  setnames(dt, 'SNP', snp_col)
  setnames(dt, snp_col, 'SNP')
  per_panel[[lbl]] <- dt
  cat(sprintf('panel %s: %d SNPs from %s (col=%s)\n', lbl, nrow(dt), fp, ec))
}

# Intersection
snp_sets <- lapply(per_panel, function(d) d$SNP)
intersection <- Reduce(intersect, snp_sets)
cat(sprintf('\nintersection size: %d SNPs\n', length(intersection)))

# Per-panel drops
summary_rows <- lapply(names(per_panel), function(lbl) {
  n0 <- length(per_panel[[lbl]]$SNP)
  ni <- length(intersection)
  data.table(panel = lbl, n_panel = n0,
             n_intersect = ni, n_dropped_from_panel = n0 - ni,
             frac_kept = ni / n0)
})
summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, 'intersect_summary.csv'))
writeLines(intersection, file.path(out_dir, 'intersect_snplist.txt'))

# Per-panel weights restricted to the intersection, aligned in SNP order.
aligned <- lapply(per_panel, function(d) {
  d <- d[SNP %in% intersection]
  setkey(d, SNP)
  d[.(intersection), on = 'SNP']  # ensure identical order
})
saveRDS(aligned, file.path(out_dir, 'intersect_weights.rds'))

# Pairwise weight-vector correlations on the intersection (M5-style; not final
# — this is a plumbing check; the actual M5 will handle allele orientation and
# use PGS scoring rather than raw weight column).
if (length(aligned) >= 2L) {
  pn <- names(aligned)
  cat('\nPairwise weight-vector correlations on intersection:\n')
  for (i in 1:(length(pn) - 1L)) {
    for (j in (i + 1L):length(pn)) {
      x <- aligned[[i]]$effect
      y <- aligned[[j]]$effect
      ok <- is.finite(x) & is.finite(y)
      rho <- if (sum(ok) >= 3L) suppressWarnings(cor(x[ok], y[ok])) else NA_real_
      cat(sprintf('  %s vs %s: cor = %+.4f  (n=%d)\n',
                  pn[i], pn[j], rho, sum(ok)))
    }
  }
}

cat('\nwrote:\n  ', file.path(out_dir, 'intersect_snplist.txt'),
    '\n  ', file.path(out_dir, 'intersect_summary.csv'),
    '\n  ', file.path(out_dir, 'intersect_weights.rds'), '\n')
