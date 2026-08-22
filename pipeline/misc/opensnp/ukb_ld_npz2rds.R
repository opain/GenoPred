#!/usr/bin/env Rscript
# Convert per-block npz produced by ukb_ld_convert.py into R RDS files
# matching Round 8's meld_ld/ block structure (fields: chr, block_id,
# block_bp, SNP, BP, cM, A1_canon, A2_canon, R_{POP}, f_{POP}, v_{POP}).

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages({
  library(reticulate)
})

OUT     <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp/meld_ld_ukb/chr22'
NPZ_DIR <- file.path(OUT, '_npz')

# Point reticulate at the magenpy env's Python (has numpy).
py_env <- '/home/claude/micromamba/envs/magenpy'
Sys.setenv(RETICULATE_PYTHON = file.path(py_env, 'bin', 'python'))
np <- import('numpy', convert = TRUE)

POPS <- c('EUR','EAS','AFR','CSA','AMR','MID')

files <- sort(list.files(NPZ_DIR, pattern = '^block_.*\\.npz$', full.names = TRUE))
cat(sprintf('found %d npz files\n', length(files)))

n_written <- 0L
for (fp in files) {
  arr <- np$load(fp, allow_pickle = TRUE)
  # arr$files lists the keys
  keys <- as.character(arr$files)

  bl <- list()
  # Scalars and small vectors
  bl$chr      <- as.integer(arr$get('chr'))
  bl$block_id <- as.integer(arr$get('block_id'))
  bl$block_bp <- as.integer(arr$get('block_bp'))
  bl$SNP      <- as.character(arr$get('SNP'))
  bl$BP       <- as.integer(arr$get('BP'))
  bl$cM       <- as.numeric(arr$get('cM'))
  bl$A1_canon <- as.character(arr$get('A1_canon'))
  bl$A2_canon <- as.character(arr$get('A2_canon'))

  # Skip empty or tiny blocks (matches m2_core.R behaviour of dropping <30)
  if (length(bl$SNP) < 2L) {
    next
  }

  for (P in POPS) {
    R <- arr$get(sprintf('R_%s', P))
    bl[[sprintf('R_%s', P)]] <- as.matrix(R)
    bl[[sprintf('f_%s', P)]] <- as.numeric(arr$get(sprintf('f_%s', P)))
    bl[[sprintf('v_%s', P)]] <- as.numeric(arr$get(sprintf('v_%s', P)))
  }

  out <- file.path(OUT, sprintf('block_%04d.rds', bl$block_id))
  saveRDS(bl, out)
  n_written <- n_written + 1L
  if (n_written %% 10L == 0L) cat(sprintf('  %d/%d\n', n_written, length(files)))
}
cat(sprintf('wrote %d RDS files to %s\n', n_written, OUT))
