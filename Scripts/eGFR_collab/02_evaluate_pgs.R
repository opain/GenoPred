#!/usr/bin/Rscript
# Evaluate PGS against a continuous phenotype, ancestry-stratified.
#
# Three analyses (each emits a TSV; combine for plotting with 03_plot.R):
#
#   default              SumStatTune: per (method, source, target), R from a
#                        regression of the pseudo-validated (sum-stats-derived)
#                        PGS, no target tuning.
#                        Output: <output>.sumstat_tune.tsv
#
#   --indiv_tune         IndivTune: per (method, source, target), R from
#                        nested k-fold CV in the target sample. Single-source
#                        configurations: inner CV picks best hyperparameter,
#                        outer fold scores it. EUR+AFR configurations:
#                        independent inner CV picks best per-source column,
#                        outer fold joint-fits both. For PRS-CSx the outer
#                        fold additionally includes META_phi_auto (a per-SNP
#                        IVW that scalar mix cannot reproduce). Requires
#                        --pgs to be the full grid (full_grid.sscore).
#                        Output: <output>.indiv_tune.tsv
#
#   --all_columns        PRSice-style sensitivity: R for every PGS column in
#                        --pgs (e.g. each P-threshold for ptclump_only).
#                        Output: <output>.all_columns.tsv
#
# Required:
#   --pgs         PLINK2 .sscore (must contain *_AVG columns)
#   --pheno       TSV: IID <pheno> (or FID IID <pheno>)
#   --covar       TSV: IID age sex PC1 ...
#   --catalogue   catalogue.tsv from the bundle
#   --target_pop  EUR or AFR
#   --output      Output prefix
# Recommended:
#   --keep        PLINK keep file (FID IID) restricting to one target ancestry.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--pgs",        type = "character"),
  make_option("--pheno",      type = "character"),
  make_option("--covar",      type = "character"),
  make_option("--catalogue",  type = "character"),
  make_option("--target_pop", type = "character"),
  make_option("--output",     type = "character"),
  make_option("--keep",       type = "character", default = NULL),
  make_option("--pheno_name", type = "character", default = NULL),
  make_option("--indiv_tune", action = "store_true", default = FALSE),
  make_option("--all_columns",action = "store_true", default = FALSE),
  make_option("--folds",      type = "integer",   default = 5L),
  make_option("--boot",       type = "integer",   default = 1000L),
  make_option("--seed",       type = "integer",   default = 42L)
)
opt <- parse_args(OptionParser(option_list = opt_list))
set.seed(opt$seed)

# ---- Load --------------------------------------------------------------
pgs   <- fread(opt$pgs)
pheno <- fread(opt$pheno)
covar <- fread(opt$covar)
cat   <- fread(opt$catalogue, na.strings = c("", "NA"))

avg_cols <- grep("_AVG$", names(pgs), value = TRUE)
if (length(avg_cols) == 0) stop("No *_AVG columns in --pgs; rerun PLINK2 --score with header output.")
pgs_cols <- sub("_AVG$", "", avg_cols)
setnames(pgs, avg_cols, pgs_cols)
pgs <- pgs[, c("#IID", pgs_cols), with = FALSE]
setnames(pgs, "#IID", "IID")

id_cols <- c("FID", "IID", "#IID", "#FID")
pheno_name <- if (!is.null(opt$pheno_name)) opt$pheno_name else setdiff(names(pheno), id_cols)[1]
stopifnot(!is.na(pheno_name), pheno_name %in% names(pheno))
message("Phenotype column: ", pheno_name)

dat <- merge(pheno[, c("IID", pheno_name), with = FALSE], covar, by = "IID")
dat <- merge(dat, pgs, by = "IID")
dat <- dat[complete.cases(dat[, c(pheno_name, setdiff(names(covar), c("FID", "IID"))), with = FALSE])]
message(sprintf("N after merge & complete cases: %d", nrow(dat)))

if (!is.null(opt$keep)) {
  keep <- fread(opt$keep, header = "auto")
  iid <- if (ncol(keep) >= 2) keep[[2]] else keep[[1]]
  n_before <- nrow(dat)
  dat <- dat[IID %in% iid]
  message(sprintf("N after --keep filter: %d (from %d)", nrow(dat), n_before))
  if (nrow(dat) == 0) stop("No samples remain after --keep filter; check ID format.")
} else {
  message("NOTE: --keep not supplied. PGS evaluation should be ancestry-stratified.")
}

covar_terms <- setdiff(names(covar), c("FID", "IID"))
y       <- dat[[pheno_name]]
X_null  <- as.matrix(dat[, covar_terms, with = FALSE])
target  <- opt$target_pop

# ---- Core R helpers -----------------------------------------------------

# Compute partial R between PGS-component (linear combination of PGS terms
# from the full model) and phenotype, controlling for covariates.
# Returns signed correlation (cor(pgs_pred, y_resid)).
fit_R <- function(yv, Xn, Xp) {
  # Xp: matrix with >=1 PGS column
  if (any(apply(Xp, 2, function(c) sd(c, na.rm = TRUE) == 0))) return(NA_real_)
  fit_full <- lm.fit(cbind(1, Xn, Xp), yv)
  fit_null <- lm.fit(cbind(1, Xn),     yv)
  pgs_part <- fit_full$fitted.values - fit_null$fitted.values
  y_resid  <- yv - fit_null$fitted.values
  if (sd(pgs_part, na.rm = TRUE) == 0) return(NA_real_)
  cor(pgs_part, y_resid)
}

sumstat_R <- function(pgs_col_names) {
  Xp <- scale(as.matrix(dat[, pgs_col_names, with = FALSE]))
  R  <- fit_R(y, X_null, Xp)
  bvals <- replicate(opt$boot, {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    fit_R(y[idx], X_null[idx, , drop = FALSE], Xp[idx, , drop = FALSE])
  })
  list(R = R, se = sd(bvals, na.rm = TRUE))
}

# Nested k-fold CV. Returns R from pooled outer-fold predictions and a
# bootstrap SE on those predictions.
pick_best_col <- function(train_dt, candidates, inner_fold_vec) {
  vals <- numeric(length(candidates)); names(vals) <- candidates
  for (col in candidates) {
    x <- scale(train_dt[[col]])[, 1]
    if (sd(x, na.rm = TRUE) == 0) { vals[col] <- -Inf; next }
    rs <- numeric(opt$folds)
    for (kk in seq_len(opt$folds)) {
      tr <- train_dt[inner_fold_vec != kk]; va <- train_dt[inner_fold_vec == kk]
      Xn_tr <- as.matrix(tr[, covar_terms, with = FALSE]); x_tr <- scale(tr[[col]])[, 1]
      ff <- lm.fit(cbind(1, Xn_tr, x_tr), tr[[pheno_name]])
      fn <- lm.fit(cbind(1, Xn_tr),       tr[[pheno_name]])
      Xn_va <- as.matrix(va[, covar_terms, with = FALSE]); x_va <- scale(va[[col]])[, 1]
      pred_full <- cbind(1, Xn_va, x_va) %*% ff$coefficients
      pred_null <- cbind(1, Xn_va)       %*% fn$coefficients
      y_va <- va[[pheno_name]]; sst <- sum((y_va - mean(y_va))^2)
      rs[kk] <- (1 - sum((y_va - pred_full)^2) / sst) -
                (1 - sum((y_va - pred_null)^2) / sst)
    }
    vals[col] <- mean(rs, na.rm = TRUE)
  }
  names(which.max(vals))
}

cv_R <- function(eur_candidates, afr_candidates = NULL, fixed_col = NULL) {
  # Single-source: pass eur_candidates only
  # Multi-source: pass both
  # fixed_col: optional column always included in the joint fit (e.g.
  #   PRS-CSx META_phi_auto, which is a per-SNP IVW of EUR/AFR betas and
  #   lives in a richer model space than a scalar EUR+AFR mix).
  fold <- sample(rep_len(seq_len(opt$folds), nrow(dat)))
  pred_full <- rep(NA_real_, nrow(dat))
  pred_null <- rep(NA_real_, nrow(dat))
  picks_eur <- character(opt$folds); picks_afr <- character(opt$folds)
  for (k in seq_len(opt$folds)) {
    tr <- dat[fold != k]; te_idx <- which(fold == k); te <- dat[fold == k]
    inner <- sample(rep_len(seq_len(opt$folds), nrow(tr)))
    p_e <- pick_best_col(tr, eur_candidates, inner); picks_eur[k] <- p_e
    p_a <- if (!is.null(afr_candidates)) pick_best_col(tr, afr_candidates, inner) else NA_character_
    picks_afr[k] <- if (is.na(p_a)) "" else p_a

    Xn_tr <- as.matrix(tr[, covar_terms, with = FALSE])
    Xp_tr <- matrix(scale(tr[[p_e]])[, 1], ncol = 1)
    if (!is.na(p_a))         Xp_tr <- cbind(Xp_tr, scale(tr[[p_a]])[, 1])
    if (!is.null(fixed_col)) Xp_tr <- cbind(Xp_tr, scale(tr[[fixed_col]])[, 1])
    ff <- lm.fit(cbind(1, Xn_tr, Xp_tr), tr[[pheno_name]])
    fn <- lm.fit(cbind(1, Xn_tr),        tr[[pheno_name]])

    Xn_te <- as.matrix(te[, covar_terms, with = FALSE])
    Xp_te <- matrix(scale(te[[p_e]])[, 1], ncol = 1)
    if (!is.na(p_a))         Xp_te <- cbind(Xp_te, scale(te[[p_a]])[, 1])
    if (!is.null(fixed_col)) Xp_te <- cbind(Xp_te, scale(te[[fixed_col]])[, 1])
    pred_full[te_idx] <- cbind(1, Xn_te, Xp_te) %*% ff$coefficients
    pred_null[te_idx] <- cbind(1, Xn_te)        %*% fn$coefficients
  }
  pgs_part <- pred_full - pred_null
  y_resid  <- y - pred_null
  R <- cor(pgs_part, y_resid, use = "complete.obs")
  bvals <- replicate(opt$boot, {
    idx <- sample.int(length(pgs_part), length(pgs_part), replace = TRUE)
    cor(pgs_part[idx], y_resid[idx], use = "complete.obs")
  })
  list(R = R, se = sd(bvals, na.rm = TRUE),
       picks_EUR = paste(unique(picks_eur), collapse = ","),
       picks_AFR = paste(unique(picks_afr[picks_afr != ""]), collapse = ","))
}

# ---- Enumerate (method, source, target) configurations ------------------
# Methods to consider for the SumStat/Indiv frames (excluding raw _multi
# rows — those provide the EUR+AFR SumStat for their underlying method).
single_source_methods <- setdiff(unique(cat$method),
                                 c(grep("_multi$", unique(cat$method), value = TRUE),
                                   "prscsx", "xwing"))

build_configs <- function(mode = c("sumstat", "indiv")) {
  mode <- match.arg(mode)
  configs <- list()
  for (m in single_source_methods) {
    # Single-source: EUR and AFR
    for (src in c("EUR", "AFR")) {
      if (mode == "sumstat") {
        cols <- cat[method == m & source_population == src & is_pseudovalidated == TRUE &
                    column_name %in% pgs_cols, column_name]
        if (length(cols) >= 1) configs[[length(configs) + 1L]] <-
          list(method = m, source = src, kind = "sumstat_single", cols = cols)
      } else {
        cols <- cat[method == m & source_population == src & column_name %in% pgs_cols, column_name]
        if (length(cols) >= 1) configs[[length(configs) + 1L]] <-
          list(method = m, source = src, kind = "indiv_single", eur = cols, afr = NULL)
      }
    }
    # EUR+AFR
    if (mode == "sumstat") {
      # LEOPARD-combined column for this target population
      leo_method <- paste0(m, "_multi")
      cols <- cat[method == leo_method & target_population == target & column_name %in% pgs_cols,
                  column_name]
      if (length(cols) == 1L) configs[[length(configs) + 1L]] <-
        list(method = m, source = "EUR+AFR", kind = "sumstat_multi", cols = cols)
    } else {
      eur <- cat[method == m & source_population == "EUR" & column_name %in% pgs_cols, column_name]
      afr <- cat[method == m & source_population == "AFR" & column_name %in% pgs_cols, column_name]
      if (length(eur) >= 1 && length(afr) >= 1) configs[[length(configs) + 1L]] <-
        list(method = m, source = "EUR+AFR", kind = "indiv_joint", eur = eur, afr = afr)
    }
  }
  # PRS-CSx (EUR+AFR only). META_phi_auto is a per-SNP IVW combination of
  # EUR/AFR betas; a scalar mix of PGS_EUR and PGS_AFR cannot reproduce it.
  # In IndivTune we therefore always include META in the joint fit alongside
  # the inner-CV-picked EUR_phi_* and AFR_phi_* columns, so IndivTune is at
  # least as expressive as SumStat (which uses META alone).
  if ("prscsx" %in% unique(cat$method)) {
    meta <- cat[method == "prscsx" & source_population == "META" &
                is_pseudovalidated == TRUE & column_name %in% pgs_cols, column_name]
    if (mode == "sumstat") {
      if (length(meta) == 1L) configs[[length(configs) + 1L]] <-
        list(method = "prscsx", source = "EUR+AFR", kind = "sumstat_multi", cols = meta)
    } else {
      eur <- cat[method == "prscsx" & source_population == "EUR" & column_name %in% pgs_cols, column_name]
      afr <- cat[method == "prscsx" & source_population == "AFR" & column_name %in% pgs_cols, column_name]
      if (length(eur) >= 1 && length(afr) >= 1) configs[[length(configs) + 1L]] <-
        list(method = "prscsx", source = "EUR+AFR", kind = "indiv_joint",
             eur = eur, afr = afr, fixed = if (length(meta) == 1L) meta else NULL)
    }
  }
  configs
}

# ---- All-columns sensitivity (PRSice style) -----------------------------
if (opt$all_columns) {
  cat_rel <- cat[(is.na(target_population) | target_population == target) &
                 column_name %in% pgs_cols]
  message(sprintf("All-columns mode: %d columns", nrow(cat_rel)))
  rows <- list()
  for (col in cat_rel$column_name) {
    res <- sumstat_R(col)
    cat_row <- cat_rel[column_name == col]
    rows[[length(rows) + 1L]] <- data.table(
      column_name = col, method = cat_row$method, source_gwas = cat_row$source_gwas,
      source_population = cat_row$source_population, hyperparameter = cat_row$hyperparameter,
      target_population = target,
      R = res$R, R_se = res$se
    )
  }
  out <- rbindlist(rows); setorder(out, -R)
  fwrite(out, paste0(opt$output, ".all_columns.tsv"), sep = "\t")
  message("Wrote: ", opt$output, ".all_columns.tsv")
  quit(save = "no")
}

# ---- SumStatTune --------------------------------------------------------
if (!opt$indiv_tune) {
  configs <- build_configs("sumstat")
  message(sprintf("SumStatTune: %d configurations", length(configs)))
  rows <- list()
  for (cfg in configs) {
    res <- sumstat_R(cfg$cols)
    rows[[length(rows) + 1L]] <- data.table(
      method = cfg$method, source = cfg$source, target = target,
      tune = "SumStatTune", R = res$R, R_se = res$se,
      pgs_cols = paste(cfg$cols, collapse = ";")
    )
  }
  out <- rbindlist(rows)
  fwrite(out, paste0(opt$output, ".sumstat_tune.tsv"), sep = "\t")
  message("Wrote: ", opt$output, ".sumstat_tune.tsv")
  quit(save = "no")
}

# ---- IndivTune ----------------------------------------------------------
configs <- build_configs("indiv")
message(sprintf("IndivTune: %d configurations, %d-fold nested CV", length(configs), opt$folds))
rows <- list()
for (cfg in configs) {
  message(sprintf("  %s / %s ...", cfg$method, cfg$source))
  res <- if (cfg$kind == "indiv_single") cv_R(cfg$eur)
         else cv_R(cfg$eur, cfg$afr, fixed_col = cfg$fixed)
  rows[[length(rows) + 1L]] <- data.table(
    method = cfg$method, source = cfg$source, target = target,
    tune = "IndivTune", R = res$R, R_se = res$se,
    picks_EUR = res$picks_EUR, picks_AFR = res$picks_AFR
  )
}
out <- rbindlist(rows)
fwrite(out, paste0(opt$output, ".indiv_tune.tsv"), sep = "\t")
message("Wrote: ", opt$output, ".indiv_tune.tsv")
