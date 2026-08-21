#!/usr/bin/env Rscript
# Round 4: build a private refdir hijacking the AMR slot with a MELD-tailored
# resample of 1KG+HGDP EUR + AFR at effective p_eur = 0.576 / p_afr = 0.424
# (the AF composition of the w = 0.50 synthetic GWAS from Round 3).
#
# Strategy:
#   1. Symlink every shared-refdir file/dir into the private refdir
#   2. OVERRIDE keep_files/AMR.keep with our resampled subset
#   3. OVERRIDE freq_files/AMR/ref.AMR.chr*.afreq recomputed on that subset
#
# Everything else (per-pop keep_files, freq_files, ref.pop.txt, per-chr genotype)
# is left symlinked so the pipeline's validation loop is happy.

suppressPackageStartupMessages(library(data.table))

# Which population slot to hijack: 'AMR' (default; uses EUR LD blocks in dbslmm) or
# 'AFR' (diagnostic; uses AFR LD blocks in dbslmm — coherent with the 42% AFR resample).
#
# Usage:
#   Rscript build_meld_refdir.R [<HIJACK_POP> [<w_eur> <n_total>]]
#
# Modes:
#   1. No w_eur      → Round 4 default: 665 EUR + 489 AFR = 1154 (p_eur = 0.576).
#                      Output dir: meld_test_r4_refdir_{hijack}
#   2. w_eur given   → Round 6 invariant-4 cross-check: draw n_eur_take = round(w_eur * n_total)
#                      and n_afr_take = round((1 - w_eur) * n_total), capped at 665 EUR
#                      and 688 AFR available.
#                      Output dir: meld_test_r6_refdir_{hijack}_w{XX}   where XX = 100*w_eur.
args <- commandArgs(trailingOnly = TRUE)
hijack_pop <- if (length(args) >= 1) args[[1]] else 'AMR'
stopifnot(hijack_pop %in% c('EUR','EAS','AFR','CSA','AMR','MID'))

if (length(args) >= 3L) {
  w_eur   <- as.numeric(args[[2]])
  n_total <- as.integer(args[[3]])
  stopifnot(is.finite(w_eur), w_eur >= 0, w_eur <= 1, n_total > 0)
  mode <- 'r6'
} else {
  w_eur   <- NA_real_
  n_total <- NA_integer_
  mode <- 'r4'
}

shared <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref'
if (mode == 'r4') {
  private <- sprintf('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r4_refdir_%s',
                     tolower(hijack_pop))
} else {
  private <- sprintf('/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r6_refdir_%s_w%02d',
                     tolower(hijack_pop), as.integer(round(w_eur * 100)))
}
plink2 <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_/bin/plink2'
seed   <- 2026L

# Available counts in the shared reference (avoid hard-coding 665/688 downstream).
avail_eur <- length(readLines(file.path(shared, 'keep_files', 'EUR.keep')))
avail_afr <- length(readLines(file.path(shared, 'keep_files', 'AFR.keep')))
stopifnot(avail_eur == 665L, avail_afr == 688L)

if (mode == 'r4') {
  n_eur_take <- 665L  # all EUR
  n_afr_take <- 489L  # gives 665/1154 = 0.576, 489/1154 = 0.424
} else {
  # Ideal counts; may exceed availability, in which case we clip and rescale
  # (the other pop's count picks up the slack in the sense that the resulting
  # composition still hits n_total exactly, but at a w slightly off nominal).
  ideal_eur <- round(w_eur       * n_total)
  ideal_afr <- round((1 - w_eur) * n_total)
  n_eur_take <- min(ideal_eur, avail_eur)
  n_afr_take <- min(ideal_afr, avail_afr)
  achieved_w <- n_eur_take / (n_eur_take + n_afr_take)
  cat(sprintf('mode=r6, requested w_eur=%.4f (n_total=%d)\n', w_eur, n_total))
  cat(sprintf('  ideal: n_eur=%d n_afr=%d\n', ideal_eur, ideal_afr))
  cat(sprintf('  taken: n_eur=%d n_afr=%d (achieved w_eur=%.4f)\n',
              n_eur_take, n_afr_take, achieved_w))
  if (ideal_eur > avail_eur || ideal_afr > avail_afr) {
    warning(sprintf('capped by availability (max EUR=%d, AFR=%d); achieved w_eur=%.4f',
                    avail_eur, avail_afr, achieved_w))
  }
}

cat('=== hijacking population slot:', hijack_pop, '===\n')
cat('=== private refdir:', private, '===\n')

if (dir.exists(private)) {
  stop('private refdir already exists: ', private,
       ' — remove before running to force fresh build')
}
dir.create(private, recursive = TRUE)
dir.create(file.path(private, 'keep_files'))
dir.create(file.path(private, 'freq_files'))
dir.create(file.path(private, 'freq_files', hijack_pop))

# ---- 1. Symlink everything from shared refdir --------------------------
cat('=== symlinking shared refdir contents ===\n')

# Top-level files/dirs (except keep_files, freq_files which we partially override)
for (f in list.files(shared, full.names = FALSE)) {
  if (f %in% c('keep_files', 'freq_files', 'genopred_1kg_hgdp.tar.gz')) next
  file.symlink(file.path(shared, f), file.path(private, f))
}

# keep_files: symlink each per-pop file except the one we hijack
for (p in setdiff(c('EUR','EAS','AFR','CSA','AMR','MID'), hijack_pop)) {
  file.symlink(file.path(shared, 'keep_files', paste0(p, '.keep')),
               file.path(private, 'keep_files', paste0(p, '.keep')))
}

# freq_files: symlink each per-pop dir except the hijacked (dir already created)
for (p in setdiff(c('EUR','EAS','AFR','CSA','AMR','MID','TRANS'), hijack_pop)) {
  file.symlink(file.path(shared, 'freq_files', p),
               file.path(private, 'freq_files', p))
}

# ---- 2. Draw the resample ----------------------------------------------
cat('\n=== drawing resample (seed =', seed, ') ===\n')
set.seed(seed)

eur_ids <- readLines(file.path(shared, 'keep_files', 'EUR.keep'))
afr_ids <- readLines(file.path(shared, 'keep_files', 'AFR.keep'))
stopifnot(length(eur_ids) == 665, length(afr_ids) == 688)

if (n_eur_take == length(eur_ids)) {
  eur_pick <- eur_ids
} else {
  eur_pick <- sample(eur_ids, n_eur_take)
}
if (n_afr_take == length(afr_ids)) {
  afr_pick <- afr_ids
} else {
  afr_pick <- sample(afr_ids, n_afr_take)
}

meld_keep <- c(eur_pick, afr_pick)
cat(sprintf('  EUR: %d (of %d), AFR: %d (of %d), total: %d\n',
            length(eur_pick), length(eur_ids),
            length(afr_pick), length(afr_ids), length(meld_keep)))
cat(sprintf('  p_eur = %.4f, p_afr = %.4f\n',
            length(eur_pick)/length(meld_keep),
            length(afr_pick)/length(meld_keep)))

writeLines(meld_keep, file.path(private, 'keep_files', paste0(hijack_pop, '.keep')))

# ---- 3. Recompute per-chr AFs on the resample --------------------------
cat(sprintf('\n=== recomputing freq_files/%s/*.afreq via plink2 ===\n', hijack_pop))

keep_path <- file.path(private, 'keep_files', paste0(hijack_pop, '.keep'))

for (i in 1:22) {
  out_prefix <- file.path(private, 'freq_files', hijack_pop,
                          sprintf('ref.%s.chr%d', hijack_pop, i))
  cmd <- sprintf(
    '%s --pfile %s/ref.chr%d --keep %s --freq --out %s > /dev/null 2>&1',
    plink2, shared, i, keep_path, out_prefix
  )
  status <- system(cmd)
  n_snps <- if (file.exists(paste0(out_prefix, '.afreq')))
    length(readLines(paste0(out_prefix, '.afreq'))) - 1L else 0L
  cat(sprintf('  chr%02d: status=%d, SNPs=%d\n', i, status, n_snps))
}

# ---- 4. Verification ---------------------------------------------------
cat('\n=== verification ===\n')
n_keep <- length(readLines(file.path(private, 'keep_files', paste0(hijack_pop, '.keep'))))
cat(sprintf('  keep_files/%s.keep: %d IIDs\n', hijack_pop, n_keep))

afreq_files <- list.files(file.path(private, 'freq_files', hijack_pop),
                          pattern = '\\.afreq$')
cat(sprintf('  freq_files/%s/*.afreq: %d files\n', hijack_pop, length(afreq_files)))

# Symlinks resolve?
res <- Sys.readlink(file.path(private, 'ref.chr1.pgen'))
cat('  symlink ref.chr1.pgen ->', res, '\n')

cat('\nDONE\n')
