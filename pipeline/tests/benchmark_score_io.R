#!/usr/bin/env Rscript
# Run from the repository root: Rscript pipeline/tests/benchmark_score_io.R
#   <score.gz> <reference-prefix-ending-in-chr> <shared-scratch-dir>
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('Expected score.gz, reference prefix, and scratch directory')
library(data.table)
source('functions/pgs.R')
ref <- read_pvar(args[2])[, c('CHR','SNP','A1','A2'), with=FALSE]
root <- tempfile('score-io-benchmark-', tmpdir=args[3]); dir.create(root)
tryCatch({
  map <- file.path(root, 'reference.tsv')
  fwrite(ref, map, sep='\t', col.names=FALSE, quote=FALSE)
  cat('Variant validation:\n')
  print(system.time(validate_score_reference(args[1], ref, record=FALSE)))
  cat('One-pass sharding:\n')
  print(system.time(prepare_score_shards(args[1], map, file.path(root, 'shard'), 1, 'functions/score_shards.awk')))
  files <- list.files(root, pattern='^shard', full.names=TRUE)
  cat('Variants:', nrow(ref), 'shards:', length(files), 'compressed bytes:', file.info(args[1])$size,
      'shard bytes:', sum(file.info(files)$size), 'reference map bytes:', file.info(map)$size, '\n')
  cat('Previous chromosome-by-chromosome extraction, including content comparison:\n')
  print(system.time(for(chr in unique(ref$CHR)) {
    index <- file.path(root, 'index.txt')
    write.table(format(which(ref$CHR == chr) + 1, scientific=FALSE), index,
                row.names=FALSE, col.names=FALSE, quote=FALSE)
    output <- file.path(root, paste0('old.', chr))
    cmd <- paste('gzip -cd --', shQuote(args[1]), '| awk',
                 shQuote('NR==FNR {rows[$1]; next} FNR==1 || FNR in rows'),
                 shQuote(index), '- | cut -d\' \' --complement -f1-3 | sed',
                 shQuote('1 s/SCORE_/score_file_1./g'), '>', shQuote(output))
    stopifnot(system2('bash', c('-o', 'pipefail', '-c', shQuote(cmd))) == 0)
    stopifnot(identical(readLines(output), readLines(file.path(root, paste0('shard.', chr, '.txt')))))
  }))
}, finally=unlink(root, recursive=TRUE))
