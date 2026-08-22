#!/bin/bash
# Round 7b Section 4 helper: for each candidate label (MELD_lambda0 + 5 single-pops),
# create a private ldpred2_ldref directory whose AMR subfolder points at the
# candidate's built LDpred2 refdir, and write a matching GenoPred config yaml.
#
# Prerequisites:
#   - build_ldpred2_ref_from_meld.R has been run for each candidate, producing
#     /users/.../pipeline/misc/opensnp/ldpred2_ref_meld/<label>/{LD_with_blocks_chr*.rds,map.rds}
#
# Output:
#   - /users/.../misc/opensnp/private_ldpred2_ldref_<label>/AMR/  → symlink to built ref
#   - /users/.../misc/opensnp/config_meld_r7bS4_ldpred2_<label>.yaml

set -euo pipefail

BASE=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/misc/opensnp
BUILT=$BASE/ldpred2_ref_meld

LABELS="MELD_lambda0 EUR EAS AFR CSA AMR equal"

for L in $LABELS; do
  PRIV=$BASE/private_ldpred2_ldref_${L}
  rm -rf "$PRIV"
  mkdir -p "$PRIV"
  # AMR slot = the built candidate ref
  ln -s "$BUILT/$L" "$PRIV/AMR"
  echo "$PRIV/AMR -> $BUILT/$L"

  # Write config
  CFG=$BASE/config_meld_r7bS4_ldpred2_${L}.yaml
  cat > "$CFG" <<EOF
# MELD Round 7b Section 4 M5 — LDpred2 with panel: ${L}
# Yengo-2022 height all-ancestry sumstats, AMR-slot hijack pattern.
outdir: /users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r7bS4_ldpred2_${L}
config_file: misc/opensnp/config_meld_r7bS4_ldpred2_${L}.yaml
gwas_list: misc/opensnp/gwas_list_r7bS4.txt
target_list: misc/opensnp/target_list.txt
refdir: /users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r4_refdir_amr
ldpred2_ldref: ${PRIV}
pgs_methods: ['ldpred2']
ldpred2_model: ['auto']
ldpred2_inference: T
testing: NA
restrict_to_target_variants: T
EOF
  echo "  config: $CFG"
done

echo ""
echo "DONE — 6 private refdirs and configs prepared."
