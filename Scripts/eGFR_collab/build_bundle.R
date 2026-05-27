#!/usr/bin/Rscript
# Build the eGFR PGS share bundle from a completed GenoPred run.
#
# Output:
#   <output_dir>/
#     catalogue.tsv
#     scores/pseudovalidated.score.gz
#     scores/full_grid.score.gz
#     scripts/01_compute_pgs.sh   (copied)
#     scripts/02_evaluate_pgs.R   (copied)
#     README.md                   (copied)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--config",      type = "character", help = "GenoPred config.yaml"),
  make_option("--genopred_dir",type = "character", default = NULL,
              help = "Path to GenoPred repo root [auto-detected from script location]"),
  make_option("--target_pops", type = "character", default = "EUR,AFR",
              help = "Comma-separated target populations for LEOPARD/xwing pseudo-val [%default]"),
  make_option("--output_dir",  type = "character", help = "Output bundle directory"),
  make_option("--full_grid",   type = "logical",   default = TRUE,
              help = "Also emit full_grid.score.gz [%default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

stopifnot(!is.null(opt$config), !is.null(opt$output_dir))

# Auto-detect GenoPred root from this script's location (Scripts/eGFR_collab/)
if (is.null(opt$genopred_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
  opt$genopred_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
}
message("Using GenoPred root: ", opt$genopred_dir)

# Source GenoPred helpers (find_pseudo, list_score_files, read_param, constants)
source(file.path(opt$genopred_dir, "functions", "constants.R"))
source(file.path(opt$genopred_dir, "functions", "misc.R"))
source(file.path(opt$genopred_dir, "functions", "pgs.R"))
source(file.path(opt$genopred_dir, "functions", "pipeline.R"))

outdir <- read_param(config = opt$config, param = "outdir", return_obj = FALSE, quiet = TRUE)
gwas_list <- read_param(config = opt$config, param = "gwas_list", quiet = TRUE)
gwas_groups <- read_param(config = opt$config, param = "gwas_groups", quiet = TRUE)
target_pops <- strsplit(opt$target_pops, ",")[[1]]

score_index <- list_score_files(opt$config, quiet = TRUE)
if (is.null(score_index) || nrow(score_index) == 0) stop("No completed score files found.")

dir.create(file.path(opt$output_dir, "scores"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "scripts"), recursive = TRUE, showWarnings = FALSE)

# Helper: source population(s) for a (method, name) row
src_pop <- function(method, name) {
  if (name %in% gwas_list$name) return(gwas_list$population[gwas_list$name == name])
  if (!is.null(gwas_groups) && name %in% gwas_groups$name) {
    g <- strsplit(gwas_groups$gwas[gwas_groups$name == name], ",")[[1]]
    return(paste(sort(unique(gwas_list$population[gwas_list$name %in% g])), collapse = "+"))
  }
  NA_character_
}

# Helper: pseudo-val column(s) for a (method, name).
# Methods needing target_pop are called once per target pop; results that are
# constant across target pops are collapsed to a single target-agnostic row.
needs_target_pop <- function(method) {
  grepl("_multi$", method) || method %in% pgs_group_methods
}

pseudo_cols <- function(method, name) {
  if (needs_target_pop(method)) {
    out <- data.table()
    for (tp in target_pops) {
      pv <- tryCatch(find_pseudo(config = opt$config, gwas = name, pgs_method = method,
                                 target_pop = tp, quiet = TRUE),
                     error = function(e) NA_character_)
      if (!is.na(pv)) out <- rbind(out, data.table(param = pv, target_population = tp))
    }
    if (nrow(out) > 0 && length(unique(out$param)) == 1L) {
      # Target-agnostic pseudo-val (e.g. prscsx META_phi_auto)
      out <- data.table(param = out$param[1], target_population = NA_character_)
    }
    return(out)
  }
  pv <- tryCatch(find_pseudo(config = opt$config, gwas = name, pgs_method = method, quiet = TRUE),
                 error = function(e) NA_character_)
  if (is.na(pv)) return(data.table())
  data.table(param = pv, target_population = NA_character_)
}

# Read each score file, build per-method tables + catalogue rows
catalogue <- list()
all_scores <- list()                    # method__name -> data.table(SNP, A1, A2, <renamed scores>)
join_key   <- NULL                      # data.table(SNP, A1, A2) — reference SNP ordering

for (i in seq_len(nrow(score_index))) {
  m <- score_index$method[i]; n <- score_index$name[i]
  fp <- file.path(outdir, "reference", "pgs_score_files", m, n, paste0("ref-", n, ".score.gz"))
  if (!file.exists(fp)) next
  dt <- fread(fp)
  score_cols <- grep("^SCORE_", names(dt), value = TRUE)
  if (length(score_cols) == 0) next

  if (is.null(join_key)) {
    join_key <- dt[, .(SNP, A1, A2)]
  } else if (nrow(dt) != nrow(join_key) || !all(dt$SNP == join_key$SNP) ||
             !all(dt$A1 == join_key$A1)) {
    stop("SNP ordering / alleles differ in ", fp, " — would need an outer join. Aborting.")
  }

  # Rename to globally unique names: <method>__<name>__<param>
  params <- sub("^SCORE_", "", score_cols)
  new_names <- paste(m, n, params, sep = "__")
  setnames(dt, score_cols, new_names)
  all_scores[[paste(m, n, sep = "__")]] <- dt[, c("SNP", "A1", "A2", new_names), with = FALSE]

  # Pseudo-val rows for this method/name
  pv <- pseudo_cols(m, n)
  pv_params <- pv$param

  for (j in seq_along(params)) {
    is_pv <- params[j] %in% pv_params
    tp <- if (is_pv && nrow(pv) > 0) pv$target_population[match(params[j], pv$param)] else NA_character_
    catalogue[[length(catalogue) + 1L]] <- data.table(
      column_name        = new_names[j],
      method             = m,
      source_gwas        = n,
      source_population  = src_pop(m, n),
      hyperparameter     = params[j],
      is_pseudovalidated = is_pv,
      target_population  = tp,
      description        = sprintf("%s PGS from %s (param=%s)%s",
                                   m, n, params[j],
                                   if (!is.na(tp)) sprintf("; tuned for target=%s", tp) else "")
    )
  }
}

catalogue_dt <- rbindlist(catalogue)
fwrite(catalogue_dt, file.path(opt$output_dir, "catalogue.tsv"), sep = "\t")
message(sprintf("catalogue.tsv: %d columns across %d methods",
                nrow(catalogue_dt), length(unique(catalogue_dt$method))))

# Build wide score tables
wide_full <- Reduce(function(a, b) cbind(a, b[, !c("SNP", "A1", "A2"), with = FALSE]),
                    all_scores, accumulate = FALSE)
fwrite(wide_full, file.path(opt$output_dir, "scores", "full_grid.score.gz"),
       sep = " ", compress = "gzip", na = "0")
message(sprintf("full_grid.score.gz: %d SNPs x %d score columns",
                nrow(wide_full), ncol(wide_full) - 3L))

pv_cols <- catalogue_dt[is_pseudovalidated == TRUE, column_name]
wide_pv <- wide_full[, c("SNP", "A1", "A2", pv_cols), with = FALSE]
fwrite(wide_pv, file.path(opt$output_dir, "scores", "pseudovalidated.score.gz"),
       sep = " ", compress = "gzip", na = "0")
message(sprintf("pseudovalidated.score.gz: %d SNPs x %d score columns",
                nrow(wide_pv), ncol(wide_pv) - 3L))

# Copy companion files
here <- file.path(opt$genopred_dir, "Scripts", "eGFR_collab")
file.copy(file.path(here, "01_compute_pgs.sh"), file.path(opt$output_dir, "scripts"), overwrite = TRUE)
file.copy(file.path(here, "02_evaluate_pgs.R"), file.path(opt$output_dir, "scripts"), overwrite = TRUE)
file.copy(file.path(here, "README.md"),         file.path(opt$output_dir),            overwrite = TRUE)

message("Bundle written to: ", opt$output_dir)
