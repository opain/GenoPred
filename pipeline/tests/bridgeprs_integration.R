#!/usr/bin/env Rscript
# Run in the pipeline analysis R environment from the repository root.
# Exercises actual target-scoring profile creation with a deterministic PLINK stub.
library(data.table)
library(yaml)
root <- normalizePath('.')
work <- tempfile('bridgeprs-integration-'); dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
dir.create(file.path(work, 'run'))
stopifnot(file.symlink(file.path(root, 'functions'), file.path(work, 'functions')))
out <- file.path(work, 'output')
ref <- file.path(work, 'reference')
dir.create(ref, recursive = TRUE)
writeLines(c('#CHROM\tPOS\tID\tREF\tALT', '22\t1\trs1\tG\tA'), file.path(ref, 'ref.chr22.pvar'))
artifact <- file.path(out, 'reference/pgs_score_files/bridgeprs/group1')
dir.create(artifact, recursive = TRUE)
con <- gzfile(file.path(artifact, 'ref-group1.score.gz'), 'wt')
writeLines(c('SNP A1 A2 SCORE_test', 'rs1 A G 0.5'), con); close(con)
writeLines(c('Param Mean SD', 'SCORE_test 0 1'), file.path(artifact, 'ref-group1-EUR.scale'))
writeLines(c('name population', 'g1 EUR', 'g2 AFR'), file.path(work, 'gwas.txt'))
writeLines(c('name gwas', 'group1 g1,g2'), file.path(work, 'groups.txt'))
writeLines(c('name path type indiv_report unrel', 'subject unused plink1 T NA'), file.path(work, 'targets.txt'))
config <- file.path(work, 'config.yaml')
write_yaml(list(outdir=out, refdir=ref, gwas_list=file.path(work, 'gwas.txt'),
                gwas_groups=file.path(work, 'groups.txt'), target_list=file.path(work, 'targets.txt'),
                score_list='NA', prepared_score_list='NA', pgs_methods=list('bridgeprs'),
                tlprs_methods='NA', leopard_methods='NA'), config)
plink <- file.path(work, 'plink2-stub')
writeLines(c('#!/bin/sh', 'while [ "$#" -gt 0 ]; do',
             '  case "$1" in --out) output="$2"; shift 2;; *) shift;; esac', 'done',
             "printf '#FID\\tIID\\tscore_file_1.test_SUM\\nF1\\tI1\\t1.25\\nF2\\tI2\\t-0.5\\n' > \"${output}.sscore\""), plink)
Sys.chmod(plink, '0755')
old <- setwd(file.path(work, 'run')); on.exit(setwd(old), add = TRUE)
script <- file.path(root, 'Scripts/target_scoring/target_scoring_pipeline.R')
status <- system2('Rscript', shQuote(c(script, '--config', config, '--name', 'subject',
                                        '--population', 'EUR', '--plink2', plink, '--test', 'chr22', '--n_cores', '1')))
if(status != 0) stop('Target scoring failed: ', status)
profile <- file.path(out, 'subject/pgs/EUR/bridgeprs/group1/subject-group1-EUR.profiles')
result <- fread(profile)
stopifnot(identical(result$FID, c('F1', 'F2')), identical(result$IID, c('I1', 'I2')),
          identical(result$group1_test, c(1.25, -0.5)))
cat('BridgePRS target profiles integration passed\n')
