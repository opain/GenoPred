#!/usr/bin/env bash
# Round 9 — download the six Pan-UKB HapMap3+ LD reference tarballs from
# Zenodo record 14614207 (Zabad et al.), verify md5, extract in place.
# Skips the 20.6 GB EUR_18m_variants.tar.gz (explicitly out of scope).
#
# ~1.8 GB total across 6 populations.

set -eu -o pipefail

DEST=/users/k1806347/oliverpainfel/Data/ukb/zenodo_14614207
LOG=$DEST/download.log

# Published md5 checksums, from Zenodo record metadata
declare -A MD5=(
  [EUR.tar.gz]=41826edf74f9cc14b3e97024119ad2e6
  [CSA.tar.gz]=38958331db0aba28edbd0eeec924d92a
  [AFR.tar.gz]=77e9c6c62ea36f88c894694e68611d99
  [EAS.tar.gz]=783ad30af1557875acc1ab6e7c32897e
  [MID.tar.gz]=79a8ece6420d2a0d579d4fef1626a953
  [AMR.tar.gz]=18c8130e167ce8cc404693524c74571e
)

# Recorded per-group N counts (matched to the prompt's table).
declare -A N=(
  [EUR]=362446
  [CSA]=8284
  [AFR]=6255
  [EAS]=2700
  [MID]=1567
  [AMR]=987
)

mkdir -p "$DEST"
cd "$DEST"

: > "$LOG"
echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) start ===" | tee -a "$LOG"

# 1. Panel-size CSV — recorded up front so it doesn't get relooked-up later.
{
  echo "pop,N"
  for P in EUR EAS AFR CSA AMR MID; do
    echo "$P,${N[$P]}"
  done
} > ukb_panel_sizes.csv
echo "wrote ukb_panel_sizes.csv" | tee -a "$LOG"

# 2. Download each with wget --continue; verify md5; extract.
fail=0
for f in EUR.tar.gz CSA.tar.gz AFR.tar.gz EAS.tar.gz MID.tar.gz AMR.tar.gz; do
  echo "--- $f ---" | tee -a "$LOG"
  if [[ ! -s "$f" ]]; then
    wget --continue --tries=5 --waitretry=10 --no-verbose \
      "https://zenodo.org/records/14614207/files/${f}?download=1" -O "$f" 2>&1 | tee -a "$LOG"
  else
    echo "already present" | tee -a "$LOG"
  fi

  got=$(md5sum "$f" | awk '{print $1}')
  want=${MD5[$f]}
  if [[ "$got" == "$want" ]]; then
    echo "md5 ok: $got" | tee -a "$LOG"
  else
    echo "MD5 MISMATCH: got=$got want=$want" | tee -a "$LOG"
    fail=$((fail+1))
    continue
  fi

  # Extract into per-pop subdir. Tarballs unpack ./chr_1..22 into CWD, so
  # every pop must extract into its own directory or later ones overwrite
  # earlier ones.
  pop=${f%.tar.gz}
  if [[ ! -d "$pop/chr_22" ]]; then
    mkdir -p "$pop"
    tar xzf "$f" -C "$pop" 2>&1 | tee -a "$LOG"
  else
    echo "already extracted at $pop/" | tee -a "$LOG"
  fi
done

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) done: failures=$fail ===" | tee -a "$LOG"
exit "$fail"
