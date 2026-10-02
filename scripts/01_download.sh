#!/bin/bash
# Download FAERS quarterly ASCII files 2012Q4 → latest into $FAERS_DIR/raw (default: ../data/raw next to this script; ~3 GB).
# FDA server is slow per connection (~0.5 MB/s), so fetch 8 quarters in parallel. Re-runnable: skips valid zips.
FAERS_DIR="${FAERS_DIR:-$(cd "$(dirname "$0")/.." && pwd)/data}"
mkdir -p "$FAERS_DIR/raw" && cd "$FAERS_DIR/raw" || exit 1
get() { f=faers_ascii_$1.zip
  [ -s $f ] && unzip -tq $f >/dev/null 2>&1 && return
  curl -sfL --retry 5 -o $f.part "https://fis.fda.gov/content/Exports/$f" && mv $f.part $f && unzip -tq $f >/dev/null 2>&1 \
    && echo "OK $f $(du -h $f | cut -f1)" || { echo "MISS $f"; rm -f $f $f.part; }; }
export -f get
for y in $(seq 2012 2026); do for q in 1 2 3 4; do [ $y -eq 2012 ] && [ $q -lt 4 ] && continue; echo ${y}q$q; done; done \
  | xargs -P 8 -I{} bash -c 'get {}'
echo DONE
