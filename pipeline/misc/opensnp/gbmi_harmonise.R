#!/usr/bin/env Rscript
# gbmi_harmonise.R — Round 8 Section 3, chr22 first.
#
# For one trait, apply the four filter stages from meld_round8_gbmi_prompt.md §3:
#   1. HapMap3 rsid intersection (default GenoPred refdir)
#   2. Per-arm rsid intersection (variant present in every single-ancestry arm)
#   3. Per-SNP N filter — retain if |N − median(N)| / median(N) ≤ 0.10 (from meta)
#   4. Allele harmonisation — align to HapMap3 REF/ALT orientation; drop
#      strand-ambiguous A/T + C/G where MAF > 0.4
#
# Emits one wide RDS per trait with harmonised meta stats plus per-arm AF, N,
# β, SE (needed downstream for weights, composition scalar, and cross-checks).
# Appends per-stage variant counts to gbmi_r8_variant_counts.csv.
#
# Usage:
#   Rscript gbmi_harmonise.R <trait>
# where <trait> is one of {Asthma, COPD, Gout, HF, IPF, Stroke, VTE}.

suppressPackageStartupMessages(library(data.table))

args  <- commandArgs(trailingOnly = TRUE)
TRAIT <- args[1]
if (is.na(TRAIT)) stop('usage: gbmi_harmonise.R <trait>')

CHR22_DIR <- '/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI/chr22'
REF_PVAR  <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref/ref.chr22.pvar'
OUT_DIR   <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'

TRAIT_ARMS <- list(
  Asthma = c('afr','amr','eas','eur','sas'),
  COPD   = c('afr','amr','eas','eur'),
  Gout   = c('afr','amr','eas','eur'),
  HF     = c('afr','amr','eas','eur'),
  IPF    = c('afr','amr','eas','eur'),
  Stroke = c('afr','amr','eas','eur'),
  VTE    = c('afr','amr','eas','eur')
)
if (!(TRAIT %in% names(TRAIT_ARMS)))
  stop(sprintf('Unknown trait %s. Choose from: %s', TRAIT,
               paste(names(TRAIT_ARMS), collapse = ', ')))
arms <- TRAIT_ARMS[[TRAIT]]

# Map GBMI arm code → GenoPred super-pop code (SAS ↔ CSA).
POP_MAP <- c(afr='AFR', amr='AMR', eas='EAS', eur='EUR', sas='CSA', mid='MID')

find_file <- function(trait, arm = NULL) {
  fn <- if (is.null(arm))
          sprintf('%s_Bothsex_inv_var_meta_GBMI_052021_nbbkgt1.chr22.tsv.gz', trait)
        else
          sprintf('%s_Bothsex_%s_inv_var_meta_GBMI_052021_nbbkgt1.chr22.tsv.gz', trait, arm)
  fp <- file.path(CHR22_DIR, fn)
  if (!file.exists(fp)) stop(sprintf('missing: %s', fp))
  fp
}

read_gbmi <- function(fp) {
  d <- fread(fp, showProgress = FALSE)
  setnames(d, '#CHR', 'chr')
  d
}

# ---------------------------------------------------------------------------
cat(sprintf('=== harmonise %s ===\n', TRAIT))

# 0. Load HapMap3 chr22 reference (SNPs only)
pvar <- fread(REF_PVAR, showProgress = FALSE)
setnames(pvar, '#CHROM', 'chr')
pvar_snp <- pvar[nchar(REF) == 1 & nchar(ALT) == 1]
# Deduplicate on rsid (defensive — expected unique but let's be safe)
pvar_snp <- pvar_snp[!duplicated(pvar_snp$ID)]
cat(sprintf('HapMap3 chr22: %d rows, %d SNP-only after dedup\n', nrow(pvar), nrow(pvar_snp)))

# 1. Read meta + arms
meta <- read_gbmi(find_file(TRAIT))
arm_data <- lapply(setNames(arms, arms), function(a) read_gbmi(find_file(TRAIT, a)))
cat(sprintf('meta: %d rows\n', nrow(meta)))
for (a in arms) cat(sprintf('arm[%s]: %d rows\n', a, nrow(arm_data[[a]])))

# Track per-stage variant counts
counts <- data.table(
  trait = TRAIT, stage = 'raw_meta', pop = 'meta', count = nrow(meta))
add_row <- function(dt, stage, pop, count)
  rbind(dt, data.table(trait = TRAIT, stage = stage, pop = pop, count = count))
for (a in arms)
  counts <- add_row(counts, 'raw_arm', unname(POP_MAP[a]), nrow(arm_data[[a]]))

# --- Filter 1: HapMap3 rsid intersection + drop rsid=='NA' + drop indels
hm3_rsids <- pvar_snp$ID
drop_and_dedup <- function(d) {
  d <- d[rsid != 'NA' & rsid %in% hm3_rsids & nchar(REF) == 1 & nchar(ALT) == 1]
  d[!duplicated(d$rsid)]
}
meta <- drop_and_dedup(meta)
for (a in arms) arm_data[[a]] <- drop_and_dedup(arm_data[[a]])
counts <- add_row(counts, 'stage1_hm3', 'meta', nrow(meta))
for (a in arms)
  counts <- add_row(counts, 'stage1_hm3', unname(POP_MAP[a]), nrow(arm_data[[a]]))
cat(sprintf('stage 1 (HapMap3): meta=%d, arms=%s\n', nrow(meta),
            paste(sapply(arm_data, nrow), collapse = ',')))

# --- Filter 2: rsid intersection across meta + every arm
common <- Reduce(intersect, c(list(meta$rsid), lapply(arm_data, function(a) a$rsid)))
meta <- meta[rsid %in% common]
for (a in arms) arm_data[[a]] <- arm_data[[a]][rsid %in% common]
counts <- add_row(counts, 'stage2_arm_intersect', 'all', length(common))
cat(sprintf('stage 2 (per-arm intersect): %d common rsids\n', length(common)))

# --- Filter 3: per-SNP N filter (from meta)
meta[, N := N_case + N_ctrl]
med_N <- median(meta$N)
keep <- abs(meta$N - med_N) / med_N <= 0.10
cat(sprintf('stage 3 (N filter): median N = %.0f; retained %d/%d (%.3f)\n',
            med_N, sum(keep), length(keep), sum(keep) / length(keep)))
meta <- meta[keep]
counts <- add_row(counts, 'stage3_N_filter_10pct', 'all', nrow(meta))

# Sensitivity: also compute counts at 5% and 20% (against the pre-filter N distribution
# so the sensitivity numbers are directly comparable to the raw stage-2 count).
for (tol in c(0.05, 0.20)) {
  n_here <- sum(abs(meta$N - med_N) / med_N <= tol)
  counts <- add_row(counts,
                    sprintf('stage3_N_filter_%dpct_sens', round(tol * 100)),
                    'all', n_here)
}

# Sync arms to the meta-retained rsids
for (a in arms) arm_data[[a]] <- arm_data[[a]][rsid %in% meta$rsid]

# --- Filter 4: allele harmonisation vs HapMap3
# Merge meta with pvar to attach HapMap3 REF/ALT.
meta_p <- merge(meta, pvar_snp[, .(rsid = ID, ref_hm3 = REF, alt_hm3 = ALT)],
                by = 'rsid', all.x = FALSE)
setkey(meta_p, rsid)

comp <- c(A = 'T', T = 'A', C = 'G', G = 'C')
meta_p[, ref_c := comp[REF]]
meta_p[, alt_c := comp[ALT]]
meta_p[, match_type := 'drop']
meta_p[REF   == ref_hm3 & ALT   == alt_hm3, match_type := 'direct']
meta_p[REF   == alt_hm3 & ALT   == ref_hm3, match_type := 'swap']
meta_p[ref_c == ref_hm3 & alt_c == alt_hm3, match_type := 'complement']
meta_p[ref_c == alt_hm3 & alt_c == ref_hm3, match_type := 'complement_swap']

# Strand-ambiguous check (A/T + C/G).
ambig <- (meta_p$REF == 'A' & meta_p$ALT == 'T') |
         (meta_p$REF == 'T' & meta_p$ALT == 'A') |
         (meta_p$REF == 'C' & meta_p$ALT == 'G') |
         (meta_p$REF == 'G' & meta_p$ALT == 'C')
maf_meta <- pmin(meta_p$all_meta_AF, 1 - meta_p$all_meta_AF)

# Drop strand-ambiguous with MAF > 0.4 (unresolvable).
# For lower-MAF ambig, assume no strand flip (i.e., treat as direct/swap).
drop_ambig <- ambig & maf_meta > 0.4
cat(sprintf('stage 4 (allele): ambig=%d, drop_high_MAF_ambig=%d, unmatched=%d\n',
            sum(ambig), sum(drop_ambig), sum(meta_p$match_type == 'drop')))
meta_p <- meta_p[!drop_ambig & match_type != 'drop']

# Flip β and AF where alleles are swapped.
flip <- meta_p$match_type %in% c('swap', 'complement_swap')
meta_p[flip, `:=`(inv_var_meta_beta = -inv_var_meta_beta,
                  all_meta_AF       = 1 - all_meta_AF)]

# Sanity: assert every arm carries the same REF/ALT as meta for the retained rsids.
for (a in arms) {
  x <- arm_data[[a]][match(meta_p$rsid, arm_data[[a]]$rsid)]
  stopifnot(!any(is.na(x$REF)))
  if (!all(x$REF == meta_p$REF & x$ALT == meta_p$ALT))
    stop(sprintf('arm %s: REF/ALT mismatch with meta on some rsids', a))
  arm_data[[a]] <- x
}

# Apply the same flip to each arm's β and AF.
for (a in arms) {
  arm_data[[a]][flip, `:=`(inv_var_meta_beta = -inv_var_meta_beta,
                           all_meta_AF       = 1 - all_meta_AF)]
}
counts <- add_row(counts, 'stage4_allele_harmonise', 'all', nrow(meta_p))
cat(sprintf('final: %d variants retained on chr22\n', nrow(meta_p)))

# ---------------------------------------------------------------------------
# 5. Build the harmonised wide table.
out <- meta_p[, .(rsid, chr, pos = POS,
                  ref = ref_hm3, alt = alt_hm3,
                  meta_beta = inv_var_meta_beta, meta_se = inv_var_meta_sebeta,
                  meta_p    = inv_var_meta_p,   meta_af = all_meta_AF,
                  meta_het_p = inv_var_het_p,
                  N = N, N_case, N_ctrl, n_bbk)]
for (a in arms) {
  P <- unname(POP_MAP[a])
  x <- arm_data[[a]]
  out[, sprintf('af_%s',   P) := x$all_meta_AF]
  out[, sprintf('beta_%s', P) := x$inv_var_meta_beta]
  out[, sprintf('se_%s',   P) := x$inv_var_meta_sebeta]
  out[, sprintf('N_%s',    P) := x$N_case + x$N_ctrl]
}

out_rds <- file.path(OUT_DIR, sprintf('gbmi_r8_%s_chr22.rds', TRAIT))
saveRDS(out, out_rds)
cat(sprintf('wrote %s (%d rows, %d cols)\n', out_rds, nrow(out), ncol(out)))

# Append counts CSV
counts_file <- file.path(OUT_DIR, 'gbmi_r8_variant_counts.csv')
if (file.exists(counts_file)) {
  prior <- fread(counts_file)
  # Overwrite this trait's rows on rerun.
  prior <- prior[trait != TRAIT]
  counts <- rbind(prior, counts)
}
fwrite(counts, counts_file)
cat(sprintf('appended to %s\n', counts_file))

cat('DONE\n')
