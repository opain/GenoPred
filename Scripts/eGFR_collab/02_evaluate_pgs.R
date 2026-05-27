#!/usr/bin/Rscript
# Evaluate PGS against a continuous phenotype.
#
# Default mode (recommended):
#   For each *pseudo-validated* PGS column relevant to the target population,
#   fit  pheno ~ covariates + scale(PGS)  and report incremental R^2 vs the
#   covariates-only model, with a bootstrap 95% CI. No hyperparameter tuning,
#   no selection bias.
#
# CV mode (--cv):
#   Nested k-fold CV across all PGS columns in --pgs to estimate an unbiased
#   "best PGS across methods" incremental R^2. Inner folds pick the best
#   column; outer folds score it on held-out data.
#
# Required inputs:
#   --pgs         PLINK2 .sscore file from 01_compute_pgs.sh
#                 (must contain *_AVG columns, one per score column)
#   --pheno       Phenotype TSV with columns: FID IID <pheno>
#   --covar       Covariate TSV with columns: FID IID age sex PC1 PC2 ...
#                 (sex coded numerically, e.g. 0/1)
#   --catalogue   catalogue.tsv from the bundle
#   --target_pop  Target population label, e.g. EUR or AFR
#                 (filters which PGS columns are relevant)
#   --output      Output TSV prefix
# Optional:
#   --pheno_name  Name of phenotype column in --pheno (default: 3rd column)
#   --cv          Run nested CV across all PGS columns
#   --folds       k for outer/inner folds (default 5)
#   --boot        Bootstrap replicates for CI (default 1000)

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
  make_option("--pheno_name", type = "character", default = NULL),
  make_option("--cv",         action = "store_true", default = FALSE),
  make_option("--folds",      type = "integer", default = 5L),
  make_option("--boot",       type = "integer", default = 1000L),
  make_option("--seed",       type = "integer", default = 42L)
)
opt <- parse_args(OptionParser(option_list = opt_list))
set.seed(opt$seed)

# ---- Load data ------------------------------------------------------------
pgs   <- fread(opt$pgs)
pheno <- fread(opt$pheno)
covar <- fread(opt$covar)
cat   <- fread(opt$catalogue, na.strings = c("", "NA"))

# PLINK2 stores per-score columns suffixed _AVG (or _SUM if cols=+scoresums)
avg_cols <- grep("_AVG$", names(pgs), value = TRUE)
if (length(avg_cols) == 0) stop("No *_AVG columns in --pgs; rerun PLINK2 --score with header output.")
# Strip _AVG suffix to recover catalogue column_name
pgs_cols <- sub("_AVG$", "", avg_cols)
setnames(pgs, avg_cols, pgs_cols)
pgs <- pgs[, c("#IID", pgs_cols), with = FALSE]
setnames(pgs, "#IID", "IID")

# Phenotype: name = user-specified, else the first non-ID column
id_cols <- c("FID", "IID", "#IID", "#FID")
pheno_name <- if (!is.null(opt$pheno_name)) {
  opt$pheno_name
} else {
  setdiff(names(pheno), id_cols)[1]
}
stopifnot(!is.na(pheno_name), pheno_name %in% names(pheno))
message("Phenotype column: ", pheno_name)

# Merge by IID (collaborator may not have FID; tolerate either)
dat <- merge(pheno[, c("IID", pheno_name), with = FALSE], covar, by = "IID")
dat <- merge(dat, pgs, by = "IID")
dat <- dat[complete.cases(dat[, c(pheno_name, setdiff(names(covar), c("FID", "IID"))), with = FALSE])]
message(sprintf("N after merge & complete cases: %d", nrow(dat)))

covar_terms <- setdiff(names(covar), c("FID", "IID"))
covar_rhs   <- paste(covar_terms, collapse = " + ")

# ---- Helpers --------------------------------------------------------------
incr_r2 <- function(y, X_full, X_null) {
  # Returns incremental R^2 = R^2(full) - R^2(null), via OLS.
  r2 <- function(X) {
    f <- lm.fit(cbind(1, as.matrix(X)), y)
    1 - sum(f$residuals^2) / sum((y - mean(y))^2)
  }
  r2(X_full) - r2(X_null)
}

boot_ci <- function(y, X_full, X_null, B = opt$boot) {
  n <- length(y)
  vals <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    vals[b] <- incr_r2(y[idx], X_full[idx, , drop = FALSE], X_null[idx, , drop = FALSE])
  }
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}

# ---- Filter catalogue to columns relevant for this target population ------
cat_rel <- cat[is_pseudovalidated == TRUE &
                 (is.na(target_population) | target_population == opt$target_pop)]
cat_rel <- cat_rel[column_name %in% pgs_cols]
message(sprintf("Pseudo-validated columns relevant for target_pop=%s: %d",
                opt$target_pop, nrow(cat_rel)))

# ---- Default mode: per-PGS incremental R^2 --------------------------------
y <- dat[[pheno_name]]
X_null <- as.matrix(dat[, covar_terms, with = FALSE])

results <- data.table()
for (col in cat_rel$column_name) {
  x <- scale(dat[[col]])[, 1]
  if (sd(x, na.rm = TRUE) == 0) next
  X_full <- cbind(X_null, PGS = x)
  r2 <- incr_r2(y, X_full, X_null)
  ci <- boot_ci(y, X_full, X_null)
  # also get per-PGS beta + p-value for reference
  m <- summary(lm(y ~ X_full))
  beta <- m$coefficients["X_fullPGS", "Estimate"]
  pval <- m$coefficients["X_fullPGS", "Pr(>|t|)"]
  results <- rbind(results, data.table(
    column_name = col, incr_R2 = r2, R2_lower = ci[1], R2_upper = ci[2],
    beta = beta, p_value = pval
  ))
}
results <- merge(results, cat_rel[, .(column_name, method, source_gwas,
                                       source_population, target_population)],
                 by = "column_name", sort = FALSE)
setorder(results, -incr_R2)
fwrite(results, paste0(opt$output, ".per_pgs.tsv"), sep = "\t")
message("Wrote: ", opt$output, ".per_pgs.tsv")

# ---- CV mode: nested k-fold across all columns ----------------------------
if (opt$cv) {
  all_cols <- intersect(pgs_cols, cat$column_name)
  message(sprintf("Nested CV across %d PGS columns, %dx%d folds",
                  length(all_cols), opt$folds, opt$folds))

  fold <- sample(rep_len(seq_len(opt$folds), nrow(dat)))
  outer_r2 <- numeric(opt$folds)
  outer_pick <- character(opt$folds)

  for (k in seq_len(opt$folds)) {
    train <- dat[fold != k]; test <- dat[fold == k]
    # Inner CV: pick best column on training data via mean inner-fold incr R^2
    inner_fold <- sample(rep_len(seq_len(opt$folds), nrow(train)))
    score_per_col <- numeric(length(all_cols)); names(score_per_col) <- all_cols
    for (col in all_cols) {
      x_tr <- scale(train[[col]])[, 1]; if (sd(x_tr, na.rm = TRUE) == 0) next
      vals <- numeric(opt$folds)
      for (kk in seq_len(opt$folds)) {
        tr2 <- train[inner_fold != kk]; va <- train[inner_fold == kk]
        Xn <- as.matrix(tr2[, covar_terms, with = FALSE])
        x2 <- scale(tr2[[col]])[, 1]
        fit_full <- lm(tr2[[pheno_name]] ~ Xn + x2)
        fit_null <- lm(tr2[[pheno_name]] ~ Xn)
        Xn_va <- as.matrix(va[, covar_terms, with = FALSE])
        x_va  <- scale(va[[col]])[, 1]
        pred_full <- cbind(1, Xn_va, x_va) %*% coef(fit_full)
        pred_null <- cbind(1, Xn_va)       %*% coef(fit_null)
        y_va <- va[[pheno_name]]
        sse_full <- sum((y_va - pred_full)^2); sse_null <- sum((y_va - pred_null)^2)
        sst <- sum((y_va - mean(y_va))^2)
        vals[kk] <- (1 - sse_full / sst) - (1 - sse_null / sst)
      }
      score_per_col[col] <- mean(vals, na.rm = TRUE)
    }
    pick <- names(which.max(score_per_col))
    outer_pick[k] <- pick

    # Outer evaluation: refit on full train, score on test
    Xn_tr <- as.matrix(train[, covar_terms, with = FALSE]); x_tr <- scale(train[[pick]])[, 1]
    fit_full <- lm(train[[pheno_name]] ~ Xn_tr + x_tr)
    fit_null <- lm(train[[pheno_name]] ~ Xn_tr)
    Xn_te <- as.matrix(test[, covar_terms, with = FALSE]); x_te <- scale(test[[pick]])[, 1]
    pred_full <- cbind(1, Xn_te, x_te) %*% coef(fit_full)
    pred_null <- cbind(1, Xn_te)       %*% coef(fit_null)
    y_te <- test[[pheno_name]]
    sse_full <- sum((y_te - pred_full)^2); sse_null <- sum((y_te - pred_null)^2)
    sst <- sum((y_te - mean(y_te))^2)
    outer_r2[k] <- (1 - sse_full / sst) - (1 - sse_null / sst)
  }
  cv_out <- data.table(
    fold = seq_len(opt$folds), picked_column = outer_pick, outer_incr_R2 = outer_r2
  )
  fwrite(cv_out, paste0(opt$output, ".cv_per_fold.tsv"), sep = "\t")
  fwrite(data.table(mean_incr_R2 = mean(outer_r2), sd = sd(outer_r2)),
         paste0(opt$output, ".cv_summary.tsv"), sep = "\t")
  message(sprintf("CV mean incremental R^2 = %.4f (sd %.4f)",
                  mean(outer_r2), sd(outer_r2)))
}
