# Shared adapters for PUMAS-EN input preparation, candidate model selection,
# and conversion back to GenoPred's canonical score-file schema.

pumas_prepare_sumstats <- function(dat) {
  raw <- as.data.frame(dat, stringsAsFactors = FALSE)
  names(raw) <- toupper(trimws(names(raw)))
  required <- c("CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "P", "N")
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("PUMAS input is missing columns: ", paste(missing, collapse = ", "))

  frequency_col <- if ("MAF" %in% names(raw)) "MAF" else if ("FREQ" %in% names(raw)) "FREQ" else if ("REF.FREQ" %in% names(raw)) "REF.FREQ" else NA_character_
  if (is.na(frequency_col)) stop("PUMAS input requires MAF, FREQ, or REF.FREQ")
  numeric_column <- function(name) suppressWarnings(as.numeric(as.character(raw[[name]])))
  frequency <- numeric_column(frequency_col)
  maf <- if (frequency_col == "MAF") frequency else pmin(frequency, 1 - frequency)

  prepared <- data.table::data.table(
    CHR = numeric_column("CHR"),
    BP = numeric_column("BP"),
    SNP = trimws(as.character(raw$SNP)),
    A1 = toupper(trimws(as.character(raw$A1))),
    A2 = toupper(trimws(as.character(raw$A2))),
    MAF = maf,
    BETA = numeric_column("BETA"),
    SE = numeric_column("SE"),
    P = numeric_column("P"),
    N = numeric_column("N")
  )

  valid <- grepl("^rs[0-9]+$", prepared$SNP, ignore.case = TRUE) &
    is.finite(prepared$CHR) & prepared$CHR %in% 1:22 &
    is.finite(prepared$BP) & prepared$BP > 0 &
    nzchar(prepared$A1) & nzchar(prepared$A2) & prepared$A1 != prepared$A2 &
    is.finite(prepared$MAF) & prepared$MAF > 0 & prepared$MAF <= 0.5 &
    is.finite(prepared$BETA) & is.finite(prepared$SE) & prepared$SE > 0 &
    is.finite(prepared$P) & prepared$P >= 0 & prepared$P <= 1 &
    is.finite(prepared$N) & prepared$N > 0
  valid[is.na(valid)] <- FALSE
  prepared <- prepared[valid]
  if (!nrow(prepared)) stop("No usable rsID summary statistics remain for PUMAS-EN")

  data.table::setorder(prepared, SNP, P)
  prepared <- unique(prepared, by = "SNP")
  data.table::setorder(prepared, CHR, BP, SNP)
  prepared[, P := pmax(P, .Machine$double.xmin)]
  prepared[, .(CHR, BP, SNP, A1, A2, MAF, BETA, SE, P, N)]
}

pumas_read_score <- function(path) {
  score <- data.table::fread(path)
  required <- c("SNP", "A1", "A2")
  missing <- setdiff(required, names(score))
  if (length(missing)) stop("Score artifact ", path, " is missing: ", paste(missing, collapse = ", "))
  score[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)))]
  if (anyDuplicated(score$SNP)) stop("Duplicate SNP IDs in score artifact: ", path)
  score
}

pumas_common_candidates <- function(score_paths) {
  if (length(score_paths) < 2) stop("PUMAS-EN candidate selection needs full-GWAS and fold artifacts")
  tables <- lapply(score_paths, pumas_read_score)
  candidates <- names(tables[[1]])[grepl("^SCORE_", names(tables[[1]]))]
  if (!length(candidates)) stop("No SCORE_* candidate columns in ", score_paths[[1]])

  retained <- character()
  dropped <- list()
  for (candidate in candidates) {
    absent <- vapply(tables, function(tab) !candidate %in% names(tab), logical(1))
    if (any(absent)) {
      dropped[[candidate]] <- "column absent from one or more fold artifacts"
      next
    }
    numeric <- vapply(tables, function(tab) is.numeric(tab[[candidate]]), logical(1))
    if (!all(numeric)) {
      dropped[[candidate]] <- "column is not numeric in every artifact"
      next
    }
    finite <- vapply(tables, function(tab) all(is.finite(tab[[candidate]])), logical(1))
    if (!all(finite)) {
      dropped[[candidate]] <- "column contains non-finite weights in one or more artifacts"
      next
    }
    retained <- c(retained, candidate)
  }
  if (!length(retained)) stop("No common finite SCORE_* candidate remains across full and fold artifacts")
  list(columns = retained, dropped = dropped)
}

pumas_write_weight_table <- function(score_path, output_path, columns) {
  score <- pumas_read_score(score_path)
  missing <- setdiff(columns, names(score))
  if (length(missing)) stop("Candidate columns disappeared from ", score_path, ": ", paste(missing, collapse = ", "))
  if (any(!vapply(score[, ..columns], is.numeric, logical(1)))) stop("Non-numeric candidate in ", score_path)
  if (any(!is.finite(as.matrix(score[, ..columns])))) stop("Non-finite candidate weight in ", score_path)
  data.table::fwrite(score[, c("SNP", "A1", columns), with = FALSE], output_path,
                     sep = "\t", quote = FALSE, col.names = FALSE)
}

# Map PUMAS EN weights onto GenoPred's reference panel, as map_score() does for
# every built-in method: one row per reference SNP (for the chromosomes in
# play), in reference order, 0 where PUMAS gave no weight, sign-flipped where
# PUMAS's effect allele is the reference's other allele. ref_scoring.R relies
# on exactly this layout - it selects rows by reference line number and splits
# columns on single spaces.
pumas_map_to_reference <- function(weights_path, ref_plink_chr, chroms) {
  weights <- data.table::fread(weights_path)
  if (!all(c("SNP", "A1", "EN") %in% names(weights))) {
    stop("PUMAS output must contain SNP, A1, and EN columns")
  }
  weights[, `:=`(SNP = as.character(SNP), A1 = toupper(as.character(A1)), EN = as.numeric(EN))]
  if (anyDuplicated(weights$SNP)) stop("Duplicate SNP IDs in PUMAS ensemble output")
  if (any(!is.finite(weights$EN))) stop("Non-finite EN weights in PUMAS output")

  # Same column mapping as GenoPred's read_pvar(): A1 = ALT, A2 = REF.
  ref <- data.table::rbindlist(lapply(chroms, function(chr) {
    pvar <- data.table::fread(paste0(ref_plink_chr, chr, ".pvar"), colClasses = "character")
    data.table::data.table(SNP = pvar[[3]], A1 = toupper(pvar[[5]]), A2 = toupper(pvar[[4]]))
  }))
  if (!nrow(ref)) stop("GenoPred reference panel is empty: ", ref_plink_chr)

  index <- match(ref$SNP, weights$SNP)
  score <- numeric(nrow(ref))
  same <- !is.na(index) & weights$A1[index] == ref$A1
  swapped <- !is.na(index) & weights$A1[index] == ref$A2
  score[same] <- weights$EN[index[same]]
  score[swapped] <- -weights$EN[index[swapped]]
  mapped <- sum(same | swapped)
  if (!mapped) stop("No PUMAS EN weights matched the GenoPred reference panel")
  message(sprintf(
    "Mapped %d of %d PUMAS EN weights onto %d reference SNPs (%d allele-incompatible, %d absent)",
    mapped, nrow(weights), nrow(ref), sum(!is.na(index)) - mapped,
    nrow(weights) - sum(weights$SNP %in% ref$SNP)
  ))
  data.table::data.table(SNP = ref$SNP, A1 = ref$A1, A2 = ref$A2, SCORE_EN = score)
}
