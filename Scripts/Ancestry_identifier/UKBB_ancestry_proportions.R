#!/usr/bin/env Rscript

# Estimate ancestry proportions for report display only, following the
# individual-level method in the bigsnpr ancestry tutorial:
# https://privefl.github.io/bigsnpr/articles/ancestry.html

library(optparse)

option_list <- list(
  make_option('--target_plink_chr', type = 'character',
              help = 'Per-chromosome PLINK2 prefix, ending in .ref.chr'),
  make_option('--ref_freq', type = 'character',
              help = 'Published UK Biobank reference allele frequencies'),
  make_option('--projection', type = 'character',
              help = 'Published UK Biobank PC projection loadings'),
  make_option('--output', type = 'character',
              help = 'Output prefix'),
  make_option('--test', type = 'character', default = 'NA',
              help = 'Pipeline test mode; enhanced ancestry requires full genome'),
  make_option('--plink2', type = 'character', default = 'plink2'),
  make_option('--threads', type = 'integer', default = 1)
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c('target_plink_chr', 'ref_freq', 'projection', 'output')
missing_options <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing_options) > 0) {
  stop('Missing required options: ', paste(missing_options, collapse = ', '))
}

library(data.table)

output_file <- paste0(opt$output, '.tsv')
metadata_file <- paste0(opt$output, '.meta.tsv')
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

empty_result <- function() {
  data.table(
    FID = character(),
    IID = character(),
    population = character(),
    proportion = numeric()
  )
}

write_unavailable <- function(reason, target_ids = data.table(FID = character(), IID = character()),
                              n_variants_matched = 0L, n_variants_used = 0L,
                              n_chromosomes = 0L) {
  fwrite(empty_result(), output_file, sep = '\t')
  if (nrow(target_ids) == 0L) {
    target_ids <- data.table(FID = NA_character_, IID = NA_character_)
  }
  metadata <- target_ids[, .(
    FID,
    IID,
    status = 'unavailable',
    reason = as.character(reason),
    n_variants_matched = as.integer(n_variants_matched),
    n_variants_used = as.integer(n_variants_used),
    n_chromosomes = as.integer(n_chromosomes),
    projection_correlation = NA_real_
  )]
  fwrite(metadata, metadata_file, sep = '\t')
  message('Enhanced ancestry unavailable: ', reason)
}

chromosomes_from_test <- function(test) {
  if (is.null(test) || is.na(test) || test == 'NA') {
    return(as.character(1:22))
  }
  sub('^chr', '', test, ignore.case = TRUE)
}

main <- function() {
  if (!is.null(opt$test) && !is.na(opt$test) && opt$test != 'NA') {
    write_unavailable(
      'Enhanced ancestry inference requires a full-genome run; chromosome-22 test mode is not sufficient.'
    )
    return(invisible(NULL))
  }

  # Test mode exits before loading the optional analysis packages. This keeps a
  # mistakenly scheduled chr22 run non-fatal even before its per-rule env has
  # been rebuilt with the enhanced-ancestry dependencies.
  library(bigsnpr)
  library(bigstatsr)
  library(Matrix)
  library(quadprog)

  chromosomes <- chromosomes_from_test(opt$test)
  first_psam <- paste0(opt$target_plink_chr, chromosomes[1], '.psam')
  if (!file.exists(first_psam)) {
    stop('Target sample file not found: ', first_psam)
  }

  psam <- fread(first_psam, header = TRUE, check.names = FALSE)
  if (ncol(psam) < 2L) {
    stop('Target PSAM does not contain FID and IID columns: ', first_psam)
  }
  target_ids <- data.table(
    FID = as.character(psam[[1]]),
    IID = as.character(psam[[2]])
  )
  # Let the top-level error handler attribute its failure row to this sample,
  # so the report can show the real reason rather than a generic fallback.
  known_target_ids <<- target_ids

  message('Reading published UK Biobank ancestry reference.')
  all_freq <- fread(opt$ref_freq, check.names = FALSE)
  projection <- fread(opt$projection, check.names = FALSE)
  if (ncol(all_freq) <= 5L || ncol(projection) <= 5L) {
    stop('Reference files must contain five variant columns plus ancestry/PC columns.')
  }

  reference_keys <- all_freq[, 1:5, with = FALSE]
  setnames(reference_keys, names(reference_keys)[1:4], c('chr', 'pos', 'a0', 'a1'))
  population_names <- names(all_freq)[-(1:5)]
  projection_names <- names(projection)[-(1:5)]
  if (nrow(all_freq) != nrow(projection)) {
    stop('Reference frequency and projection files have different numbers of variants.')
  }
  if (length(projection_names) != 16L) {
    stop('Expected 16 PC projection columns, found ', length(projection_names), '.')
  }
  correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074,
                  1.099, 1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

  # The tutorial merges reference groups that are too close to distinguish
  # robustly. Preserve its labels so the report does not overstate precision.
  display_groups <- population_names
  display_groups[display_groups %in% c('Scandinavia', 'United Kingdom', 'Ireland')] <-
    'Europe (North West)'
  display_groups[display_groups %in% c('Europe (South East)', 'Europe (North East)')] <-
    'Europe (East)'
  display_levels <- unique(display_groups)

  projected_sum <- NULL
  reference_projection_sum <- NULL
  total_matched <- 0L
  total_used <- 0L
  chromosomes_used <- character()
  temp_dir <- tempfile('ukbb_ancestry_')
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  for (chr in chromosomes) {
    target_prefix <- paste0(opt$target_plink_chr, chr)
    pgen_file <- paste0(target_prefix, '.pgen')
    pvar_file <- paste0(target_prefix, '.pvar')
    if (!file.exists(pgen_file) || !file.exists(pvar_file)) {
      stop('Missing formatted target chromosome files for chromosome ', chr, '.')
    }

    message('Preparing chromosome ', chr, '.')
    bed_prefix <- file.path(temp_dir, paste0('target_chr', chr))
    status <- system2(
      opt$plink2,
      c('--pfile', target_prefix, '--make-bed', '--out', bed_prefix,
        '--threads', as.character(opt$threads)),
      stdout = TRUE,
      stderr = TRUE
    )
    if (!file.exists(paste0(bed_prefix, '.bed'))) {
      stop('PLINK2 could not convert chromosome ', chr, ' to BED format: ',
           paste(status, collapse = '\n'))
    }

    bim <- fread(paste0(bed_prefix, '.bim'), header = FALSE)
    if (ncol(bim) < 6L) {
      stop('Invalid PLINK BIM file for chromosome ', chr, '.')
    }
    target_keys <- bim[, .(
      chr = as.integer(V1),
      pos = as.integer(V4),
      a1 = toupper(as.character(V5)),
      a0 = toupper(as.character(V6)),
      beta = 1
    )]

    # snp_match returns a data.frame; convert so the column filters below are
    # evaluated inside the table (otherwise `beta` resolves to base::beta).
    matched <- as.data.table(snp_match(target_keys, reference_keys))
    if (nrow(matched) == 0L) {
      message('No UK Biobank reference variants matched on chromosome ', chr, '.')
      next
    }
    matched <- matched[!is.na(beta) & !is.na(`_NUM_ID_.ss`) & !is.na(`_NUM_ID_`)]
    total_matched <- total_matched + nrow(matched)
    if (nrow(matched) == 0L) next

    rds <- snp_readBed2(
      paste0(bed_prefix, '.bed'),
      ind.col = matched$`_NUM_ID_.ss`
    )
    obj <- snp_attach(rds)
    genotype <- obj$genotypes
    missing_by_variant <- big_counts(genotype)[4, ]
    usable <- which(missing_by_variant < 5)
    if (length(usable) == 0L) {
      message('All matched variants were too sparse on chromosome ', chr, '.')
      next
    }

    genotype_imputed <- snp_fastImputeSimple(genotype)
    pc_projection <- as.matrix(
      projection[matched$`_NUM_ID_`, -(1:5), with = FALSE]
    )[usable, , drop = FALSE]
    reference_freq <- as.matrix(
      all_freq[matched$`_NUM_ID_`, -(1:5), with = FALSE]
    )[usable, , drop = FALSE]
    good <- complete.cases(pc_projection) & complete.cases(reference_freq)
    usable <- usable[good]
    pc_projection <- pc_projection[good, , drop = FALSE]
    reference_freq <- reference_freq[good, , drop = FALSE]
    if (length(usable) == 0L) next

    beta <- matched$beta[usable]
    projected <- big_prodMat(
      genotype_imputed,
      sweep(pc_projection, 2, correction / 2, '*'),
      ind.col = usable,
      center = 1 - beta,
      scale = beta
    )
    projected <- as.matrix(projected)
    if (is.null(dim(projected))) {
      projected <- matrix(projected, nrow = nrow(target_ids))
    }
    reference_projection <- crossprod(pc_projection, reference_freq)

    if (is.null(projected_sum)) {
      projected_sum <- matrix(0, nrow = nrow(projected), ncol = ncol(projected))
      reference_projection_sum <- matrix(
        0, nrow = nrow(reference_projection), ncol = ncol(reference_projection)
      )
    }
    projected_sum <- projected_sum + projected
    reference_projection_sum <- reference_projection_sum + reference_projection
    total_used <- total_used + length(usable)
    chromosomes_used <- c(chromosomes_used, chr)
  }

  if (is.null(projected_sum) || total_used < 100L) {
    write_unavailable(
      paste0('Fewer than 100 usable variants matched the UK Biobank reference (',
             total_used, ' usable variants).'),
      target_ids,
      total_matched,
      total_used,
      length(chromosomes_used)
    )
    return(invisible(NULL))
  }

  # Solve one constrained quadratic program per target individual. X contains
  # the reference population centres in the same PC space as the target.
  x <- reference_projection_sum
  dmat <- nearPD(crossprod(x), base.matrix = TRUE)$mat
  amat <- cbind(1, diag(ncol(x)))
  bvec <- c(1, rep(0, ncol(x)))
  solution <- matrix(NA_real_, nrow = nrow(projected_sum), ncol = ncol(x))
  correlations <- rep(NA_real_, nrow(projected_sum))

  for (i in seq_len(nrow(projected_sum))) {
    y <- as.numeric(projected_sum[i, ])
    fit <- solve.QP(
      Dmat = dmat,
      dvec = as.numeric(crossprod(y, x)),
      Amat = amat,
      bvec = bvec,
      meq = 1
    )$solution
    solution[i, ] <- pmax(fit, 0)
    solution[i, ] <- solution[i, ] / sum(solution[i, ])
    predicted <- as.numeric(x %*% solution[i, ])
    if (sd(y) > 0 && sd(predicted) > 0) {
      correlations[i] <- cor(y, predicted)
    }
  }

  # Collapse the closely related reference groups after solving the QP.
  collapsed <- matrix(0, nrow = nrow(solution), ncol = length(display_levels))
  colnames(collapsed) <- display_levels
  for (j in seq_along(display_groups)) {
    collapsed[, display_groups[j]] <- collapsed[, display_groups[j]] + solution[, j]
  }

  result <- rbindlist(lapply(seq_len(nrow(collapsed)), function(i) {
    data.table(
      FID = target_ids$FID[i],
      IID = target_ids$IID[i],
      population = display_levels,
      proportion = as.numeric(collapsed[i, ])
    )
  }))
  metadata <- target_ids[, .(
    FID,
    IID,
    status = 'ok',
    reason = NA_character_,
    n_variants_matched = as.integer(total_matched),
    n_variants_used = as.integer(total_used),
    n_chromosomes = as.integer(length(chromosomes_used)),
    projection_correlation = correlations
  )]

  fwrite(result, output_file, sep = '\t')
  fwrite(metadata, metadata_file, sep = '\t')
  message('Enhanced ancestry completed with ', total_used, ' usable variants across ',
          length(chromosomes_used), ' chromosomes.')
}

known_target_ids <- data.table(FID = character(), IID = character())

tryCatch(
  main(),
  error = function(e) {
    write_unavailable(conditionMessage(e), known_target_ids)
    quit(save = 'no', status = 0)
  }
)
