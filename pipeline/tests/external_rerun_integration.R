#!/usr/bin/env Rscript
# Run in the pipeline analysis R environment from the repository root.
library(data.table)
root <- normalizePath('.')
work <- tempfile('external-rerun-'); dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
dir.create(file.path(work, 'run'))
stopifnot(file.symlink(file.path(root, 'functions'), file.path(work, 'functions')))
ref <- file.path(work, 'reference.chr')
saveRDS(data.table(CHR=22, SNP='rs123', BP_GRCh37=100, A1='A', A2='C', IUPAC='M'), paste0(ref, '22.rds'))
writeLines(c('#CHROM\tPOS\tID\tREF\tALT', '22\t100\trs123\tC\tA'), paste0(ref, '22.pvar'))
score <- file.path(work, 'input.txt')
output <- file.path(work, 'results/ref-score')
script <- file.path(root, 'Scripts/external_score_processor/external_score_processor.R')
old <- setwd(file.path(work, 'run')); on.exit(setwd(old), add = TRUE)
run <- function() system2('Rscript', shQuote(c(script, '--ref_plink_chr', ref, '--score', score, '--output', output, '--test', 'chr22')))
header <- 'rsID effect_allele other_allele effect_weight'
writeLines(c(header, 'rs123 A C 0.5'), score)
stopifnot(run() == 0, file.exists(paste0(output, '.score.gz')))
shared_log <- file.path(work, 'shared.log')
stopifnot(file.link(paste0(output, '.log'), shared_log))
shared_contents <- readLines(shared_log)
writeLines('stale', paste0(output, '.score.gz.reference.md5'))
writeLines(c(header, 'rs123 A C 0.5', 'rs999 A C 0.2', 'rs998 A C 0.2', 'rs997 A C 0.2'), score)
stopifnot(run() == 0, !file.exists(paste0(output, '.score.gz')),
          !file.exists(paste0(output, '.score.gz.reference.md5')),
          identical(readLines(shared_log), shared_contents),
          any(grepl('^Skipping', readLines(paste0(output, '.log')))))
cat('External-score rerun integration passed\n')
