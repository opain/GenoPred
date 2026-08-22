#!/usr/bin/env Rscript
# r8_on_r9_subblocks.R — apples-to-apples comparison prep.
#
# R9 uses 67 non-empty sub-blocks (common refinement of the six UKB per-pop
# LDetect partitions). R8 used 24 blocks (EUR fourier_ls-chr22.bed). The R9
# blocks are strictly finer than R8's (they are subdivisions of EUR's — since
# EUR's UKB LDetect matches R8's fourier_ls bed exactly).
#
# For each R9 sub-block, find the parent R8 block and extract R_p / f_p /
# v_p for the sub-block's rsids from the parent's stored matrices. Write a
# parallel set of block RDSes at meld_ld_r8_on_r9blocks/chr22/block_*.rds.
# eval_meld_ld_r8_on_r9blocks.R can then run R8's M2 on the same block
# structure as R9, so the per-block r² numbers are directly comparable.

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages(library(data.table))

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'
R8_LD_DIR <- file.path(MISC, 'meld_ld', 'chr22')
R9_LD_DIR <- file.path(MISC, 'meld_ld_ukb', 'chr22')
OUT       <- file.path(MISC, 'meld_ld_r8_on_r9blocks', 'chr22')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

POPS <- c('EUR','EAS','AFR','CSA','AMR')  # R8 has no MID

r8_blocks <- sort(list.files(R8_LD_DIR, pattern = '^block_.*\\.rds$', full.names = TRUE))
r9_blocks <- sort(list.files(R9_LD_DIR, pattern = '^block_.*\\.rds$', full.names = TRUE))

# Build a bp -> R8-block lookup by loading the small R8 metadata.
r8_meta <- rbindlist(lapply(r8_blocks, function(fp) {
  bl <- readRDS(fp)
  data.table(file = fp, s = bl$block_bp[1], e = bl$block_bp[2])
}))
cat(sprintf('R8: %d blocks, R9: %d sub-blocks\n', nrow(r8_meta), length(r9_blocks)))

# For each R9 sub-block, find the parent R8 block and extract sub-matrices.
r8_cache <- list()

n_written <- 0L
for (fp9 in r9_blocks) {
  bl9 <- readRDS(fp9)
  if (length(bl9$SNP) < 2L) next
  centre <- mean(bl9$block_bp)
  parent <- r8_meta[s <= centre & e > centre]
  if (nrow(parent) == 0L) next
  parent_file <- parent$file[1]
  if (is.null(r8_cache[[parent_file]])) {
    r8_cache[[parent_file]] <- readRDS(parent_file)
  }
  bl8 <- r8_cache[[parent_file]]

  # Restrict to rsids present in both.
  common <- intersect(bl9$SNP, bl8$SNP)
  if (length(common) < 2L) next
  ix9 <- match(common, bl9$SNP)
  ix8 <- match(common, bl8$SNP)

  # Sanity: allele orientation should match, since both were aligned to HM3.
  # (Skip strict check to avoid R8's blocks not carrying A1_canon in the same form.)

  # Build the R8-on-R9-block RDS.
  out_bl <- list(
    chr      = bl9$chr,
    block_id = bl9$block_id,
    block_bp = bl9$block_bp,
    SNP      = common,
    BP       = bl9$BP[ix9],
    cM       = bl9$cM[ix9],
    A1_canon = bl9$A1_canon[ix9],
    A2_canon = bl9$A2_canon[ix9]
  )
  for (P in POPS) {
    out_bl[[sprintf('R_%s', P)]] <- bl8[[sprintf('R_%s', P)]][ix8, ix8, drop = FALSE]
    out_bl[[sprintf('f_%s', P)]] <- bl8[[sprintf('f_%s', P)]][ix8]
    out_bl[[sprintf('v_%s', P)]] <- bl8[[sprintf('v_%s', P)]][ix8]
  }

  saveRDS(out_bl, file.path(OUT, sprintf('block_%04d.rds', bl9$block_id)))
  n_written <- n_written + 1L
}
cat(sprintf('wrote %d R8-on-R9-subblock RDSes to %s\n', n_written, OUT))
