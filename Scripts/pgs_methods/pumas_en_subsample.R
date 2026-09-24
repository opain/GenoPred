#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(
  make_option("--sumstats", type = "character"),
  make_option("--gwas", type = "character"),
  make_option("--work-dir", dest = "work_dir", type = "character"),
  make_option("--code-dir", dest = "code_dir", type = "character"),
  make_option("--ld-blocks", dest = "ld_blocks", type = "character"),
  make_option("--helpers", type = "character"),
  make_option("--folds", type = "integer", default = 4L),
  make_option("--partitions", type = "character", default = "0.6,0.2,0.1,0.1"),
  make_option("--threads", type = "integer", default = 4L)
)
opt <- parse_args(OptionParser(option_list = option_list))
source(opt$helpers)

if (length(opt$gwas) != 1 || !grepl("^[A-Za-z0-9_.]+$", opt$gwas)) stop("Invalid PUMAS trait name")
input_dir <- file.path(opt$work_dir, "input")
upstream_dir <- file.path(opt$work_dir, "upstream")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(upstream_dir, recursive = TRUE, showWarnings = FALSE)

cleaned <- pumas_prepare_sumstats(fread(opt$sumstats))
fwrite(cleaned, file.path(input_dir, paste0(opt$gwas, ".txt")), sep = "\t", quote = FALSE)

pumas_script <- file.path(opt$code_dir, "PUMA-ensemble.subsampling.R")
setwd(opt$code_dir)
args <- c(
  pumas_script,
  "--k", as.character(opt$folds),
  "--partitions", opt$partitions,
  "--trait_name", opt$gwas,
  "--ensemble", "EN",
  "--gwas_path", input_dir,
  "--ld_path", opt$ld_blocks,
  "--output_path", paste0(upstream_dir, .Platform$file.sep),
  "--parallel",
  "--threads", as.character(opt$threads)
)
status <- system2("Rscript", vapply(args, shQuote, character(1)))
if (!identical(status, 0L)) stop("Pinned PUMAS subsampling exited with status ", status)

required <- c(
  unlist(lapply(seq_len(opt$folds), function(fold) c(
    file.path(upstream_dir, sprintf("%s.gwas.omnibus.ite%d.txt", opt$gwas, fold)),
    file.path(upstream_dir, sprintf("%s.xty.omnibus.ite%d.txt", opt$gwas, fold))
  ))),
  file.path(upstream_dir, paste0(opt$gwas, ".omnibus.forEVAL.txt"))
)
missing <- required[!file.exists(required) | file.info(required)$size <= 0]
if (length(missing)) stop("PUMAS subsampling did not create expected files: ", paste(missing, collapse = ", "))
