#!/usr/bin/env Rscript
# Run in the pipeline analysis R environment from the repository root.
library(data.table)
library(yaml)
root <- normalizePath('.')
work <- tempfile('continuous-rerun-'); dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)
dir.create(file.path(work, 'run'))
stopifnot(file.symlink(file.path(root, 'functions'), file.path(work, 'functions')))
out <- file.path(work, 'output'); ref <- file.path(work, 'reference'); dir.create(ref)
writeLines(c('#CHROM\tPOS\tID\tREF\tALT', '22\t1\trs1\tG\tA'), file.path(ref, 'ref.chr22.pvar'))
fwrite(data.table(FID=paste0('F',1:8), IID=paste0('I',1:8), POP='EUR'), file.path(ref, 'ref.pop.txt'), sep=' ')
artifact <- file.path(out, 'reference/pgs_score_files/bridgeprs/group1'); dir.create(artifact, recursive=TRUE)
con <- gzfile(file.path(artifact, 'ref-group1.score.gz'), 'wt')
writeLines(c('SNP A1 A2 SCORE_test', 'rs1 A G 0.5'), con); close(con)
pcdir <- file.path(out, 'reference/pc_score_files/TRANS'); dir.create(pcdir, recursive=TRUE)
fwrite(data.table(FID=paste0('F',1:8), IID=paste0('I',1:8), PC1=1:8), file.path(pcdir, 'ref-TRANS-pcs.profiles'), sep=' ')
writeLines('PC1 0 1', file.path(pcdir, 'ref-TRANS-pcs.EUR.scale'))
writeLines(c('name population', 'g1 EUR', 'g2 AFR'), file.path(work, 'gwas.txt'))
writeLines(c('name gwas', 'group1 g1,g2'), file.path(work, 'groups.txt'))
config <- file.path(work, 'config.yaml')
write_yaml(list(outdir=out, refdir=ref, gwas_list=file.path(work, 'gwas.txt'),
                gwas_groups=file.path(work, 'groups.txt'), score_list='NA', prepared_score_list='NA',
                pgs_methods=list('bridgeprs'), tlprs_methods='NA', leopard_methods='NA'), config)
plink <- file.path(work, 'plink2-stub')
writeLines(c('#!/bin/sh', 'while [ "$#" -gt 0 ]; do',
             '  case "$1" in --out) output="$2"; shift 2;; *) shift;; esac', 'done',
             "printf '#FID\\tIID\\tscore_file_1.test_SUM\\nF1\\tI1\\t1\\nF2\\tI2\\t5\\nF3\\tI3\\t2\\nF4\\tI4\\t9\\nF5\\tI5\\t4\\nF6\\tI6\\t3\\nF7\\tI7\\t8\\nF8\\tI8\\t11\\n' > \"${output}.sscore\""), plink)
Sys.chmod(plink, '0755')
old <- setwd(file.path(work, 'run')); on.exit(setwd(old), add = TRUE)
script <- file.path(root, 'Scripts/ref_scoring/ref_scoring.R')
run <- function(continuous) system2('Rscript', shQuote(c(script, '--config', config, '--continuous', continuous,
                                                          '--plink2', plink, '--test', 'chr22', '--n_cores', '1')))
model <- file.path(artifact, 'ref-group1-TRANS.model.rds')
stopifnot(run('F') == 0, file.exists(file.path(artifact, 'ref-group1-EUR.profiles')), !file.exists(model))
stopifnot(run('T') == 0, file.exists(model), file.exists(file.path(artifact, 'ref-group1-TRANS.scale')),
          file.exists(file.path(artifact, 'ref-group1-TRANS.profiles')))
cat('Discrete-to-continuous reference scoring integration passed\n')
