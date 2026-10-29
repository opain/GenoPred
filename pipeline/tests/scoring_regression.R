#!/usr/bin/env Rscript
# Run from the repository root: Rscript pipeline/tests/scoring_regression.R
library(data.table)
library(parallel)
source('functions/pgs.R')
getDoParWorkers <- function() 1L
expect_error <- function(expr) {
  if(!inherits(tryCatch(force(expr), error = identity), 'error')) stop('Expected an error')
}

pcs <- data.table(FID = c('f','f','g'), IID = c('1','2','3'), PC1 = c(1,2,3))
scores <- data.table(FID = c('g','f','f'), IID = c('3','1','2'), SCORE = c(9,4,8))
mean_model <- lm(y ~ PC1, data = data.frame(y = c(2,4,6), PC1 = 1:3))
var_model <- lm(y ~ PC1, data = data.frame(y = c(0,0,0), PC1 = 1:3))
models <- list(SCORE = list(mean_model = mean_model, var_model = var_model))
a <- score_adjust(copy(scores), copy(pcs), models)
b <- score_adjust(copy(scores), copy(pcs[c(3,1,2)]), models)
stopifnot(identical(a, b), identical(a$SCORE, c(3,2,4)))
expect_error(score_adjust(copy(scores), pcs[-1], models))
expect_error(score_adjust(copy(scores), rbind(pcs, pcs[1]), models))
expect_error(check_sample_order(scores, scores[c(2,1,3)]))
expect_error(check_sample_order(scores, scores[-1]))
check_sample_order(scores, copy(scores))

ref <- data.table(SNP = c('rs1','rs2'), A1 = c('A','C'), A2 = c('G','T'))
weights <- data.table(SNP = c('rs2','rs1'), A1 = c('T','A'), A2 = c('C','G'), SCORE = c(3,2))
mapped <- map_score(ref, weights)
stopifnot(identical(mapped$SCORE, c(2,-3)))
bad <- copy(weights); bad$A2[1] <- 'A'
expect_error(map_score(ref, bad))
expect_error(map_score(ref, rbind(weights, weights[1])))
strand <- copy(weights); strand[1, `:=`(A1 = 'A', A2 = 'G')]
stopifnot(identical(map_score(ref, strand)$SCORE, c(2, -3)))
stats <- copy(weights); stats[, `:=`(BETA = c(3,2), FREQ = c(.2,.4), SCORE = NULL)]
matched <- allele_match(stats, ref)
stopifnot(matched[SNP == 'rs2', BETA] == -3, matched[SNP == 'rs2', FREQ] == .8)
stats[1, A2 := 'A']
expect_error(allele_match(stats, ref))

path <- tempfile(fileext = '.gz')
con <- gzfile(path, 'wt')
writeLines(c('SNP A1 A2 SCORE_test', 'rs1 A G 2', 'rs2 C T -3'), con)
close(con)
validate_score_reference(path, ref)
shard_dir <- tempfile('score shards '); dir.create(shard_dir)
ref_map <- file.path(shard_dir, 'reference map.tsv')
fwrite(data.table(CHR = c(1,2), ref), ref_map, sep = '\t', col.names = FALSE)
shard_prefix <- file.path(shard_dir, 'score shard')
prepare_score_shards(path, ref_map, shard_prefix, 7, 'functions/score_shards.awk')
stopifnot(identical(readLines(paste0(shard_prefix, '.1.txt')), c('score_file_7.test', '2')),
          identical(readLines(paste0(shard_prefix, '.2.txt')), c('score_file_7.test', '-3')))
fwrite(data.table(CHR = c(1,2), ref[2:1]), ref_map, sep = '\t', col.names = FALSE)
expect_error(prepare_score_shards(path, ref_map, shard_prefix, 7, 'functions/score_shards.awk'))
unlink(shard_dir, recursive = TRUE)
sidecar <- paste0(path, '.reference.md5')
stopifnot(file.exists(sidecar), identical(readLines(sidecar), reference_variant_fingerprint(ref)))
expect_error(validate_score_reference(path, ref[2:1]))
expect_error(validate_score_reference(path, ref[1]))
writeLines('wrong fingerprint', sidecar)
expect_error(validate_score_reference(path, ref))
unlink(c(path, sidecar))
prefix <- tempfile()
input <- paste0(prefix, '.score.gz'); file.create(input)
file.create(paste0(prefix, '-EUR.profiles'))
file.create(paste0(prefix, '-EUR.scale'))
stopifnot(!reference_score_needs_update(input, prefix), reference_score_needs_update(input, prefix, TRUE))
trans <- paste0(prefix, c('-TRANS.model.rds', '-TRANS.scale', '-TRANS.profiles'))
file.create(trans)
stopifnot(!reference_score_needs_update(input, prefix, TRUE))
unlink(trans[1]); stopifnot(reference_score_needs_update(input, prefix, TRUE))
unlink(c(input, paste0(prefix, c('-EUR.profiles', '-EUR.scale')), trans))

# Discovery must include BridgePRS group output, and must not revive an old
# external file when the current harmonisation log says it was skipped.
outdir <- tempfile(); dir.create(outdir)
on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
pgs_group_methods <- c('prscsx', 'xwing', 'bridgeprs')
pgs_methods_noneur <- c('ptclump')
read_param <- function(config, param, ...) switch(param,
  gwas_list = data.frame(name = c('g1','g2'), population = c('EUR','AFR')),
  pgs_methods = c('bridgeprs'),
  score_list = data.frame(name = 'external1'),
  prepared_score_list = NULL,
  gwas_groups = data.frame(name = 'group1', gwas = 'g1,g2'),
  tlprs_methods = NA, leopard_methods = NA,
  outdir = outdir)
bridge_dir <- file.path(outdir, 'reference/pgs_score_files/bridgeprs/group1')
external_dir <- file.path(outdir, 'reference/pgs_score_files/external/external1')
dir.create(bridge_dir, recursive = TRUE); dir.create(external_dir, recursive = TRUE)
file.create(file.path(bridge_dir, 'ref-group1.score.gz'), file.path(external_dir, 'ref-external1.score.gz'))
log <- file.path(external_dir, 'ref-external1.log')
writeLines(c('Skipping', 'Analysis finished at now'), log)
found <- list_score_files('unused', quiet = TRUE)
stopifnot(nrow(found) == 1, found$method == 'bridgeprs', found$name == 'group1')
writeLines('Analysis finished at now', log)
found <- list_score_files('unused', quiet = TRUE)
stopifnot(nrow(found) == 2, 'external' %in% found$method)
unlink(outdir, recursive = TRUE)
cat('Scoring regression checks passed\n')
