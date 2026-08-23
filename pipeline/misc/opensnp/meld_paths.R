# meld_paths.R — canonical locations of MELD data, LD panels, and results.
#
# The scripts under pipeline/misc/opensnp/ are code that lives in the git
# repo; their inputs and outputs are large or numerous, and live outside
# the repo under /users/k1806347/oliverpainfel/Analyses/MELD/.
#
# All MELD analysis and plotting scripts source this file to pick up the
# canonical paths, so relocating the analysis tree only requires
# editing MELD_ROOT below.

MELD_ROOT    <- '/users/k1806347/oliverpainfel/Analyses/MELD'
MELD_DATA    <- file.path(MELD_ROOT, 'data')
MELD_LD      <- file.path(MELD_ROOT, 'ld')
MELD_RESULTS <- file.path(MELD_ROOT, 'results')

# Return the results directory for a given round; auto-create if missing.
# Round labels are r4/r5/r6/r6b/r7/r7b/r8/r9/r10.
r_results <- function(round) {
  d <- file.path(MELD_RESULTS, round)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Convenience: harmonised sumstats RDSes live under data/gbmi_r8_harmonised/.
r8_harmonised <- function(trait) {
  file.path(MELD_DATA, 'gbmi_r8_harmonised',
            sprintf('gbmi_r8_%s_chr22.rds', trait))
}

# LD directory subtree names — one per subtree that scripts refer to.
meld_ld_1kg_hgdp <- function() file.path(MELD_LD, '1kg_hgdp')
meld_ld_ukb      <- function() file.path(MELD_LD, 'ukb')
meld_ld_r8_on_r9blocks <- function() file.path(MELD_LD, '1kg_hgdp_on_r9blocks')
ref_empirical    <- function() file.path(MELD_LD, 'ref_empirical')
