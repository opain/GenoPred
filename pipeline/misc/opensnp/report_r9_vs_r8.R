#!/usr/bin/env Rscript
# report_r9_vs_r8.R — Round 9 Section 5 assembler.
#
# Assembles a side-by-side table of Round 8 (1KG+HGDP) vs Round 9 (UKB) M2
# numbers on the SAME (R9) block structure, per trait, per candidate, at
# δ = 0.01. The R8 side uses gbmi_r8b_m2_ci.csv (R8 LD data evaluated on
# the R9 sub-block structure — apples to apples). The original R8 CIs from
# the coarser 24-block structure are also carried through in a separate
# column for reference.
#
# Outputs:
#   gbmi_r9_vs_r8.csv          — long-form table (trait × candidate × source)
#   gbmi_r9_vs_r8_summary.csv  — headline paired diffs per trait

.libPaths(c('/home/claude/Rlibs', .libPaths()))
suppressPackageStartupMessages(library(data.table))

MISC <- '/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp'

ci_r9  <- fread(file.path(MISC, 'gbmi_r9_m2_ci.csv'))
ci_r8b <- fread(file.path(MISC, 'gbmi_r8b_m2_ci.csv'))
ci_r8_original <- fread(file.path(MISC, 'gbmi_r8_m2_ci.csv'))

paired_r9  <- fread(file.path(MISC, 'gbmi_r9_m2_paired_ci.csv'))
paired_r8b <- fread(file.path(MISC, 'gbmi_r8b_m2_paired_ci.csv'))

# Filter to δ=0.01 headline
ci_r9   <- ci_r9  [delta == 0.01]
ci_r8b  <- ci_r8b [delta == 0.01]
ci_r8_o <- ci_r8_original[delta == 0.01]

# Tag each source
ci_r9  [, source := 'R9_UKB']
ci_r8b [, source := 'R8_1KG_on_R9blocks']
ci_r8_o[, source := 'R8_1KG_original_blocks']

long <- rbind(ci_r9, ci_r8b, ci_r8_o, use.names = TRUE)
fwrite(long, file.path(MISC, 'gbmi_r9_vs_r8.csv'))

# Headline paired diffs per trait: MELD-λ0 − EUR, MELD-λ0 − equal (composition),
# equal − EUR (blending), MELD-λ1 − MELD-λ0
#
# The R9 side uses R9 UKB numbers; the "R8" side for headline uses R8b (R8 LD
# on R9 blocks) so the comparison is like-for-like block-wise.

pair_map <- c(
  'MELD_lambda0 − EUR'       = 'MELD_lambda0_afproj − EUR',
  'MELD_lambda0 − equal (composition)' = 'MELD_lambda0_afproj − MELD_lambda0_equal',
  'equal − EUR (blending)'   = 'MELD_lambda0_equal − EUR',
  'MELD_lambda1 − MELD_lambda0' = 'MELD_lambda1_afproj − MELD_lambda0_afproj'
)

extract_pair <- function(dt, key_label, delta_val = 0.01) {
  cmp <- pair_map[key_label]
  x <- dt[delta == delta_val & comparison == cmp,
          .(trait, quantity = key_label, mean_diff, lo95, hi95)]
  x
}

r9_headline <- rbindlist(lapply(names(pair_map), extract_pair, dt = paired_r9))
r9_headline[, source := 'R9_UKB']
r8b_headline <- rbindlist(lapply(names(pair_map), extract_pair, dt = paired_r8b))
r8b_headline[, source := 'R8_1KG']

hl <- rbind(r9_headline, r8b_headline)

# Long → wide, per quantity: two columns "R8_1KG" and "R9_UKB" plus deltas
hl[, txt := sprintf('%+0.4f [%+0.4f, %+0.4f]', mean_diff, lo95, hi95)]

# For each trait × quantity, side-by-side + change
wide <- dcast(hl, trait + quantity ~ source, value.var = 'mean_diff')
wide_txt <- dcast(hl, trait + quantity ~ source, value.var = 'txt')

# Set factor ordering for readable output
hl$quantity <- factor(hl$quantity, levels = names(pair_map))
wide$quantity <- factor(wide$quantity, levels = names(pair_map))

# The R8 vs R9 difference-of-differences is NOT bootstrappable across LD
# sources; the plan explicitly asks for side-by-side reporting.
wide[, change_R9_minus_R8 := R9_UKB - R8_1KG]

setcolorder(wide, c('trait','quantity','R8_1KG','R9_UKB','change_R9_minus_R8'))
setorder(wide, trait, quantity)
fwrite(wide, file.path(MISC, 'gbmi_r9_vs_r8_summary.csv'))
cat('=== Section 5 summary (Δ = R9_UKB − R8_1KG on same sub-block structure) ===\n')
print(wide)

# Per-candidate table — R9 vs R8b M2 mean per trait × candidate
ci_r9_wide  <- dcast(ci_r9,  trait ~ candidate, value.var = 'mean')
ci_r8b_wide <- dcast(ci_r8b, trait ~ candidate, value.var = 'mean')

setnames(ci_r9_wide, setdiff(names(ci_r9_wide), 'trait'),
         paste0(setdiff(names(ci_r9_wide), 'trait'), '_R9'))
setnames(ci_r8b_wide, setdiff(names(ci_r8b_wide), 'trait'),
         paste0(setdiff(names(ci_r8b_wide), 'trait'), '_R8'))
both <- merge(ci_r8b_wide, ci_r9_wide, by = 'trait')
fwrite(both, file.path(MISC, 'gbmi_r9_vs_r8_candidates.csv'))
cat('\n=== Per-candidate M2 means (R8_1KG vs R9_UKB, on R9 sub-blocks) ===\n')
# print a subset: EUR, MELD_lambda0_afproj, equal
sel_cols <- c('trait',
              'EUR_R8', 'EUR_R9',
              'MELD_lambda0_afproj_R8', 'MELD_lambda0_afproj_R9',
              'MELD_lambda0_equal_R8', 'MELD_lambda0_equal_R9')
print(both[, ..sel_cols])

cat('\nwrote gbmi_r9_vs_r8.csv, gbmi_r9_vs_r8_summary.csv, gbmi_r9_vs_r8_candidates.csv\n')
