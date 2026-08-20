#!/bin/bash
# Round 5: build a MELD-tailored quickprs LD reference for the 665-EUR + 489-AFR
# resample (identical to Round 4). Reproduces docs/prep_quickprs_ref.Rmd:59-117 for the
# hijacked AMR slot only.
#
# Output: 7 runtime files at $OUT_DIR/
#   AMR.cors.{bim,bin,noise,root}
#   AMR.bld.ldak.quickprs.{matrix,tagging}
#   highld.snps

set -euo pipefail

CEPH=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/.snakemake/conda/3f88447533fd10040edfdcea8db853f7_
PLINK=$CEPH/bin/plink
PLINK2=$CEPH/bin/plink2
LDAK=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/software/ldak5.2/ldak5.2.linux

SHARED_REFDIR=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ref
LDAK_MAP=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_map/genetic_map_b37
LDAK_BLD=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_bld
LDAK_HIGHLD=/users/k1806347/oliverpainfel/Software/MyGit/GenoPred/pipeline/resources/data/ldak_highld/highld.txt

KEEP=/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r4_refdir_amr/keep_files/AMR.keep
OUT_DIR=/users/k1806347/oliverpainfel/Data/OpenSNP/GenoPred/meld_test_r5_quickprs_ldref/AMR
TMP=/home/claude/meld_run_logs/r5_ldak_tmp
NCORES=5

mkdir -p "$OUT_DIR"
rm -rf "$TMP"; mkdir -p "$TMP"

echo "=== 1. Per-chr bfile via plink2 --keep ==="
for i in $(seq 1 22); do
  $PLINK2 --pfile "$SHARED_REFDIR/ref.chr${i}" \
          --keep "$KEEP" \
          --make-bed \
          --out "$TMP/chr${i}" >/dev/null 2>&1
  echo "  chr${i} done: $(wc -l < $TMP/chr${i}.bim) SNPs, $(wc -l < $TMP/chr${i}.fam) samples"
done

echo ""
echo "=== 2. Merge chromosomes into ref_merge (plink1) ==="
for i in $(seq 2 22); do echo "$TMP/chr${i}"; done > "$TMP/merge_list.txt"
$PLINK --bfile "$TMP/chr1" \
       --merge-list "$TMP/merge_list.txt" \
       --make-bed \
       --out "$TMP/ref_merge" > "$TMP/merge.log" 2>&1 || {
  echo "merge failed — see $TMP/merge.log"; tail -20 "$TMP/merge.log"; exit 1;
}
echo "  ref_merge: $(wc -l < $TMP/ref_merge.bim) SNPs, $(wc -l < $TMP/ref_merge.fam) samples"

# Free per-chr bfiles
rm -f "$TMP"/chr*.{bed,bim,fam,log}

echo ""
echo "=== 3. Rewrite bim IDs to CHR:BP + insert genetic distances ==="
awk '{$2=$1":"$4; print $0}' "$TMP/ref_merge.bim" > "$TMP/tmp.bim"
mv "$TMP/tmp.bim" "$TMP/ref_merge.bim"

$PLINK --bfile "$TMP/ref_merge" \
       --cm-map "$LDAK_MAP/genetic_map_chr@_combined_b37.txt" \
       --make-bed \
       --out "$TMP/map" > "$TMP/cm.log" 2>&1
awk '{print $2, $3}' "$TMP/map.bim" > "$TMP/map.all"
awk '(NR==FNR){arr[$1]=$2;next}{print $1, $2, arr[$2], $4, $5, $6}' \
    "$TMP/map.all" "$TMP/ref_merge.bim" > "$TMP/tmp.bim"
mv "$TMP/tmp.bim" "$TMP/ref_merge.bim"
rm -f "$TMP"/map*
echo "  bim rewritten with genetic distances"

echo ""
echo "=== 4. LDAK: cut-weights ==="
$LDAK --cut-weights "$TMP/sections" \
      --bfile "$TMP/ref_merge" \
      --max-threads $NCORES > "$TMP/cut_weights.log" 2>&1
tail -3 "$TMP/cut_weights.log"

echo ""
echo "=== 5. LDAK: calc-weights-all ==="
$LDAK --calc-weights-all "$TMP/sections" \
      --bfile "$TMP/ref_merge" \
      --max-threads $NCORES > "$TMP/calc_weights.log" 2>&1
tail -3 "$TMP/calc_weights.log"

echo ""
echo "=== 6. Assemble bld/ (copy BLD annotations + weights.short as bld65) ==="
mkdir -p "$TMP/bld"
cp "$LDAK_BLD"/bld* "$TMP/bld/"
mv "$TMP/sections/weights.short" "$TMP/bld/bld65"
ls "$TMP/bld" | wc -l | xargs -I{} echo "  bld/ has {} files"

echo ""
echo "=== 7. LDAK: calc-tagging (this is slow) ==="
date '+  start=%H:%M:%S'
$LDAK --calc-tagging "$TMP/bld.ldak" \
      --bfile "$TMP/ref_merge" \
      --ignore-weights YES \
      --power -.25 \
      --annotation-number 65 \
      --annotation-prefix "$TMP/bld/bld" \
      --window-cm 1 \
      --save-matrix YES \
      --max-threads $NCORES > "$TMP/calc_tagging.log" 2>&1
date '+  end  =%H:%M:%S'
tail -3 "$TMP/calc_tagging.log"

echo ""
echo "=== 8. LDAK: calc-cors (also slow) ==="
date '+  start=%H:%M:%S'
$LDAK --calc-cors "$TMP/tmp" \
      --bfile "$TMP/ref_merge" \
      --window-cm 3 \
      --max-threads $NCORES > "$TMP/calc_cors.log" 2>&1
date '+  end  =%H:%M:%S'
tail -3 "$TMP/calc_cors.log"

echo ""
echo "=== 9. LDAK: cut-genes (high-LD regions) ==="
$LDAK --cut-genes "$TMP/highld" \
      --bfile "$TMP/ref_merge" \
      --genefile "$LDAK_HIGHLD" \
      --max-threads $NCORES > "$TMP/cut_genes.log" 2>&1
tail -3 "$TMP/cut_genes.log"

echo ""
echo "=== 10. Copy runtime files into $OUT_DIR ==="
for i in bim bin noise root; do
  cp "$TMP/tmp.cors.$i" "$OUT_DIR/AMR.cors.$i"
done
for i in matrix tagging; do
  cp "$TMP/bld.ldak.$i" "$OUT_DIR/AMR.bld.ldak.quickprs.$i"
done
cp "$TMP/highld/genes.predictors.used" "$OUT_DIR/highld.snps"
ls -lah "$OUT_DIR"

echo ""
echo "DONE"
