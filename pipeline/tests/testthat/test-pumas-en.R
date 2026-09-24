library(data.table)
library(testthat)

helper_candidates <- c(
  "Scripts/pgs_methods/pumas_en_helpers.R",
  "../Scripts/pgs_methods/pumas_en_helpers.R",
  "../../Scripts/pgs_methods/pumas_en_helpers.R",
  "../../../Scripts/pgs_methods/pumas_en_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not locate pumas_en_helpers.R")
source(helper_path)

test_that("PUMAS summary-statistic preparation derives MAF and removes invalid/duplicate rows", {
  raw <- data.table(
    CHR = c(1, 1, 1, 2), BP = c(100, 200, 300, 100),
    SNP = c("rs1", "rs1", "not_an_rsid", "rs2"),
    A1 = c("A", "A", "G", "T"), A2 = c("G", "G", "A", "C"),
    BETA = c(0.1, 0.2, 0.3, -0.2), SE = c(0.02, 0.03, 0.04, 0.05),
    P = c(0.01, 0.001, 0.02, 0), N = c(10000, 10000, 10000, 9000),
    FREQ = c(0.2, 0.8, 0.4, 0.7)
  )

  prepared <- pumas_prepare_sumstats(raw)

  expect_named(prepared, c("CHR", "BP", "SNP", "A1", "A2", "MAF", "BETA", "SE", "P", "N"))
  expect_equal(prepared$SNP, c("rs1", "rs2"))
  expect_equal(prepared$BETA[prepared$SNP == "rs1"], 0.2)
  expect_equal(prepared$MAF, c(0.2, 0.3))
  expect_gt(prepared$P[prepared$SNP == "rs2"], 0)
})

test_that("PUMAS retains only candidate columns shared and finite across all artifacts", {
  temp <- tempfile("pumas-candidates-")
  dir.create(temp)
  full <- file.path(temp, "full.tsv")
  fold1 <- file.path(temp, "fold1.tsv")
  fold2 <- file.path(temp, "fold2.tsv")
  fwrite(data.table(SNP = c("rs1", "rs2"), A1 = c("A", "C"), A2 = c("G", "T"),
                    SCORE_keep = c(0.1, 0.2), SCORE_drop = c(1, 2)), full, sep = "\t")
  fwrite(data.table(SNP = c("rs1", "rs2"), A1 = c("A", "C"), A2 = c("G", "T"),
                    SCORE_keep = c(0.3, 0.4), SCORE_drop = c(3, 4)), fold1, sep = "\t")
  fwrite(data.table(SNP = c("rs1", "rs2"), A1 = c("A", "C"), A2 = c("G", "T"),
                    SCORE_keep = c(0.5, 0.6)), fold2, sep = "\t")

  candidates <- pumas_common_candidates(c(full, fold1, fold2))

  expect_equal(candidates$columns, "SCORE_keep")
  expect_match(candidates$dropped$SCORE_drop, "absent")

  weights <- file.path(temp, "weights.txt")
  pumas_write_weight_table(full, weights, candidates$columns)
  headerless <- fread(weights, header = FALSE)
  expect_equal(names(headerless), c("V1", "V2", "V3"))
  expect_equal(headerless$V1, c("rs1", "rs2"))
  expect_equal(headerless$V2, c("A", "C"))
  expect_equal(headerless$V3, c(0.1, 0.2))
})

test_that("PUMAS EN weights map onto the GenoPred reference panel", {
  temp <- tempfile("pumas-map-")
  dir.create(temp)
  prefix <- file.path(temp, "ref.chr")
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT",
               "22\t100\trs1\tG\tA",
               "22\t200\trs2\tC\tT",
               "22\t300\trs3\tA\tG",
               "22\t400\trs4\tA\tC"), paste0(prefix, "22.pvar"))
  weights <- file.path(temp, "ensemble.txt")
  # rs1 matches ALT, rs2 is the swapped allele, rs4 is allele-incompatible,
  # rs9 is absent from the reference and rs3 has no PUMAS weight.
  fwrite(data.table(SNP = c("rs1", "rs2", "rs4", "rs9"), A1 = c("A", "C", "G", "A"),
                    EN = c(0.25, 0.5, 1, 2)), weights, sep = "\t")

  canonical <- pumas_map_to_reference(weights, prefix, 22)

  expect_named(canonical, c("SNP", "A1", "A2", "SCORE_EN"))
  expect_equal(canonical$SNP, c("rs1", "rs2", "rs3", "rs4"))
  expect_equal(canonical$A1, c("A", "T", "G", "C"))
  expect_equal(canonical$A2, c("G", "C", "A", "A"))
  expect_equal(canonical$SCORE_EN, c(0.25, -0.5, 0, 0))
})
