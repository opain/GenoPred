#!/usr/bin/env Rscript
# Parse GBMI manifest to build the Round 8 download list.
# Selection rules (from meld_round8_gbmi_prompt.md §1):
#   sex == "Bothsex"
#   Note == "all biobanks"
#   phenotype_short ∈ {Asthma, COPD, Gout, HF, IPF, Stroke, VTE}
#   For each trait: one multi-ancestry meta row + every single-ancestry arm row.
# Multi-ancestry meta rows have `ancestry` as a comma-separated list (e.g. "eas,afr,eur,amr,sas").
# Single-ancestry arm rows have `ancestry` in {afr, amr, eas, eur, sas, mid}.

suppressPackageStartupMessages(library(data.table))

MANIFEST <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI/manifest_GBMI_summary_statistics.csv'
OUT_DIR  <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI'
TARGETS  <- c('Asthma','COPD','Gout','HF','IPF','Stroke','VTE')

# Row 1 is the main header; row 2 is a sub-header describing the wget columns.
# fread with skip=1 puts the sub-header as column names; use header row 1 explicitly.
hdr <- readLines(MANIFEST, n = 2)
cat('Header row 1:\n', hdr[1], '\n', sep = '')
cat('Header row 2:\n', hdr[2], '\n', sep = '')

# Read with row 1 as header, skip row 2 (sub-header) after.
raw <- fread(MANIFEST, skip = 2, header = FALSE)
# Column names from row 1:
setnames(raw, c('phenotype','phenotype_short','sex','ancestry','Note','biobank',
                'wget_gz','wget_tbi','wget_qq','wget_manhattan'))

cat(sprintf('\nTotal data rows: %d\n', nrow(raw)))
cat('Unique phenotype_short:\n'); print(sort(unique(raw$phenotype_short)))
cat('Unique sex:\n'); print(sort(unique(raw$sex)))
cat('Unique Note:\n'); print(sort(unique(raw$Note)))

# Filter to Bothsex + all biobanks + target traits
sel <- raw[sex == 'Bothsex' & Note == 'all biobanks' & phenotype_short %in% TARGETS]
cat(sprintf('\nAfter Bothsex + all biobanks + target-trait filter: %d rows\n', nrow(sel)))

# Split into meta (comma-separated ancestry) and arm (single-ancestry) rows
sel[, is_meta := grepl(',', ancestry, fixed = TRUE)]
cat('\nMeta vs arm counts:\n'); print(table(sel$is_meta))

cat('\n=== Meta rows (one per trait) ===\n')
print(sel[is_meta == TRUE, .(phenotype_short, ancestry)])

cat('\n=== Arm rows ===\n')
print(sel[is_meta == FALSE, .(phenotype_short, ancestry)])

cat('\n=== Sanity: arms per trait ===\n')
print(sel[is_meta == FALSE, .N, by = phenotype_short])

# Extract the URL out of each wget command.
# Pattern: `wget <URL>  -O <fname>`
extract_url <- function(cmd) {
  m <- regmatches(cmd, regexpr('https?://[^ ]+', cmd))
  if (length(m) == 0) NA_character_ else m
}
sel[, url_gz  := sapply(wget_gz,  extract_url)]
sel[, url_tbi := sapply(wget_tbi, extract_url)]
sel[, fname_gz  := basename(url_gz)]
sel[, fname_tbi := basename(url_tbi)]

# Write two things:
#   (a) a tidy CSV of the 36 rows with url+fname
#   (b) a plain-text URL list for wget --input-file
out_csv <- file.path(OUT_DIR, 'gbmi_r8_download_manifest.csv')
fwrite(sel[, .(phenotype_short, ancestry, is_meta, biobank,
               url_gz, fname_gz, url_tbi, fname_tbi)], out_csv)
cat(sprintf('\nWrote %s (%d rows)\n', out_csv, nrow(sel)))

out_urls <- file.path(OUT_DIR, 'gbmi_r8_urls.txt')
writeLines(c(sel$url_gz, sel$url_tbi), out_urls)
cat(sprintf('Wrote %s (%d URLs)\n', out_urls, 2 * nrow(sel)))

# Final summary table
cat('\n=== Final download list — grouped by trait ===\n')
print(sel[order(phenotype_short, is_meta, ancestry),
          .(trait = phenotype_short, kind = ifelse(is_meta, 'META', 'ARM'),
            ancestry, fname_gz)])
