#!/usr/bin/env bash
# Round 8 GBMI downloader.
# Reads URLs from gbmi_r8_urls.txt (produced by gbmi_parse_manifest.R).
# Verifies via gzip -t / .tbi presence; logs to gbmi_r8_download.log.
# Resumable via wget -c.

set -eu -o pipefail

DEST=/users/k1806347/oliverpainfel/Data/GWAS_sumstats/GBMI
URLS=$DEST/gbmi_r8_urls.txt
LOG=$DEST/gbmi_r8_download.log

cd "$DEST"

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) start ===" | tee -a "$LOG"
echo "URL count: $(wc -l < "$URLS")" | tee -a "$LOG"

# wget: -c continue, -nc skip if already present would clash with -c; use -c only.
# --tries=5 --waitretry=10 for transient S3 glitches.
# One serial wget invocation is fine here — S3 has plenty of bandwidth per connection
# and we want a predictable log ordering. If we ever want parallelism, GNU parallel.
wget --continue --tries=5 --waitretry=10 --no-verbose \
     --input-file="$URLS" \
     -a "$LOG"

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) wget done ===" | tee -a "$LOG"

# Post-download verification
fail=0
while read -r url; do
  fn=$(basename "$url")
  if [[ ! -s "$fn" ]]; then
    echo "MISSING or empty: $fn" | tee -a "$LOG"
    fail=$((fail+1))
    continue
  fi
  # gzip integrity check only on .gz (not on .tbi which is bgzip-format but tabix will validate)
  case "$fn" in
    *.gz)
      if ! gzip -t "$fn" 2>>"$LOG"; then
        echo "GZIP CORRUPT: $fn" | tee -a "$LOG"
        fail=$((fail+1))
      fi
      ;;
  esac
done < "$URLS"

echo "=== $(date -u +%Y-%m-%dT%H:%M:%SZ) verification done, failures=$fail ===" | tee -a "$LOG"
exit "$fail"
