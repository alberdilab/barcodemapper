#!/usr/bin/env bash
set -euo pipefail
source "$BENCH_ROOT/tasks/1.2.3-comparators/common.sh"
OUT="$RES/minimap2-sensitive"
if [ "$TOOL" = vsearch ]; then FMT=vsearch; EXT=tsv; else FMT=paf; EXT=paf; fi
for f in accepted.$EXT hitmap.tsv calls.tsv manifest.json SHA256SUMS; do
  [ ! -e "$OUT/$RUN.$f" ] || { echo "exists, not overwriting: $OUT/$RUN.$f" >&2; exit 4; }
done
: > "$OUT/$RUN.accepted.$EXT"
for m in 1 2; do
  python -I "$TASKS/accept_hits.py" $FMT $m "$OUT/$RUN.R$m.raw.$EXT" "$OUT/$RUN.R$m.accepted.$EXT" | tee "$OUT/$RUN.R$m.accept.log"
  cat "$OUT/$RUN.R$m.accepted.$EXT" >> "$OUT/$RUN.accepted.$EXT"
done
zcat "$REF" | grep '^>' | python -I "$TASKS/make_hit_map.py" $FMT "$OUT/$RUN.accepted.$EXT" "$OUT/$RUN.hitmap.tsv"
python -I "$TASKS/convert_calls.py" $TOOL --read-ids "$READS/${RUN}_1.fq.gz" \
  --input "$OUT/$RUN.accepted.$EXT" --accession-paths "$OUT/$RUN.hitmap.tsv" --output "$OUT/$RUN.calls.tsv"
python -I "$TASKS/write_manifest_sensitive.py" "$BENCH_ROOT" "$RUN" "$TOOL"
