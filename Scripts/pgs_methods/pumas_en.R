#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--gwas", type = "character"),
  make_option("--methods", type = "character"),
  make_option("--full-scores", dest = "full_scores", type = "character"),
  make_option("--fold-scores", dest = "fold_scores", type = "character"),
  make_option("--folds", type = "integer", default = 4L),
  make_option("--partitions", type = "character", default = "0.6,0.2,0.1,0.1"),
  make_option("--threads", type = "integer", default = 4L),
  make_option("--code-dir", dest = "code_dir", type = "character"),
  make_option("--work-dir", dest = "work_dir", type = "character"),
  make_option("--reference-prefix", dest = "reference_prefix", type = "character"),
  make_option("--ref-plink-chr", dest = "ref_plink_chr", type = "character"),
  make_option("--test", type = "character", default = "NA"),
  make_option("--out-score", dest = "out_score", type = "character"),
  make_option("--helpers", type = "character"),
  make_option("--pumas-commit", dest = "pumas_commit", type = "character"),
  make_option("--resource-manifest", dest = "resource_manifest", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
source(opt$helpers)

split_csv <- function(value) trimws(strsplit(value, ",", fixed = TRUE)[[1]])
methods <- split_csv(opt$methods)
full_scores <- split_csv(opt$full_scores)
fold_scores <- split_csv(opt$fold_scores)
if (length(methods) < 2 || length(full_scores) != length(methods) ||
    length(fold_scores) != length(methods) * opt$folds) {
  stop("PUMAS method, full-score, and fold-score lists have inconsistent lengths")
}
if (anyDuplicated(methods)) stop("Duplicate component methods passed to PUMAS-EN")
trait <- opt$gwas
work_dir <- opt$work_dir
fold_weights <- file.path(work_dir, "weights", "folds")
full_weights <- file.path(work_dir, "weights", "full")
xty_dir <- file.path(work_dir, "upstream")
stats_dir <- file.path(work_dir, "upstream")
eval_dir <- file.path(work_dir, "evaluation")
for (dir in c(fold_weights, full_weights, eval_dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

candidate_info <- list()
for (method_index in seq_along(methods)) {
  method <- methods[[method_index]]
  fold_paths <- vapply(seq_len(opt$folds), function(fold) {
    fold_scores[[((fold - 1L) * length(methods)) + method_index]]
  }, character(1))
  paths <- c(full_scores[[method_index]], fold_paths)
  candidates <- pumas_common_candidates(paths)
  candidate_info[[method]] <- list(columns = candidates$columns, dropped = candidates$dropped)

  pumas_write_weight_table(
    full_scores[[method_index]],
    file.path(full_weights, paste0(trait, ".", method, ".ite1.txt")),
    candidates$columns
  )
  for (fold in seq_len(opt$folds)) {
    pumas_write_weight_table(
      fold_paths[[fold]],
      file.path(fold_weights, paste0(trait, ".", method, ".ite", fold, ".txt")),
      candidates$columns
    )
  }
}

pumas_script <- file.path(opt$code_dir, "PUMAS-ensemble.evaluation.R")
eval_args <- c(
  pumas_script,
  "--k", as.character(opt$folds),
  "--ref_path", opt$reference_prefix,
  "--trait_name", trait,
  "--prs_method", paste(methods, collapse = ","),
  "--ensemble", "EN",
  "--xty_path", paste0(xty_dir, .Platform$file.sep),
  "--stats_path", paste0(stats_dir, .Platform$file.sep),
  "--weight_path", paste0(fold_weights, .Platform$file.sep),
  "--full_weight_path", paste0(full_weights, .Platform$file.sep),
  "--output_path", paste0(eval_dir, .Platform$file.sep),
  "--parallel",
  "--threads", as.character(opt$threads)
)
setwd(opt$code_dir)
status <- system2("Rscript", vapply(eval_args, shQuote, character(1)))
if (!identical(status, 0L)) stop("Pinned PUMAS evaluation exited with status ", status)

weights_path <- file.path(eval_dir, paste0(trait, ".ensemble.weights.txt"))
if (!file.exists(weights_path) || file.info(weights_path)$size <= 0) stop("PUMAS evaluation did not write ensemble weights")
chroms <- if (is.na(opt$test) || opt$test == "NA") 1:22 else as.numeric(gsub("chr", "", opt$test))
canonical <- pumas_map_to_reference(weights_path, opt$ref_plink_chr, chroms)
dir.create(dirname(opt$out_score), recursive = TRUE, showWarnings = FALSE)
fwrite(canonical, opt$out_score, sep = " ", quote = FALSE, compress = "gzip")

r2_files <- list.files(eval_dir, pattern = "\\.en\\.r2\\.txt$", full.names = TRUE)
r2 <- lapply(r2_files, function(path) scan(path, what = numeric(), quiet = TRUE))
names(r2) <- basename(r2_files)
provenance <- list(
  method = "PUMAS-EN",
  ensemble = "EN",
  pumas_commit = opt$pumas_commit,
  pumas_reference_sha256 = readLines(opt$resource_manifest, warn = FALSE),
  gwas = trait,
  components = methods,
  candidate_columns = lapply(candidate_info, `[[`, "columns"),
  dropped_candidates = lapply(candidate_info, `[[`, "dropped"),
  folds = opt$folds,
  partitions = as.numeric(split_csv(opt$partitions)),
  per_fold_en_r2 = r2
)
write_json(provenance, file.path(dirname(opt$out_score), "pumas_en.provenance.json"),
           auto_unbox = TRUE, pretty = TRUE, null = "null")
