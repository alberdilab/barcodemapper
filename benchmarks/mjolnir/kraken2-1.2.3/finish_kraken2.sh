#!/usr/bin/env bash
set -euo pipefail
for SETTING in kraken2 kraken2-conf0.1; do
  for RUN in SRR14351851 ERR6001161 ERR6001173; do
    export SETTING RUN
    source "$BENCH_ROOT/tasks/1.2.3-kraken2/common.sh"
    for f in calls.tsv manifest.json SHA256SUMS; do
      [ ! -e "$OUT/$RUN.$f" ] || { echo "exists, not overwriting: $OUT/$RUN.$f" >&2; exit 4; }
    done
    python -I "$TASKS/convert_calls.py" kraken2 --read-ids "$READS/${RUN}_1.fq.gz" \
      --input "$OUT/$RUN.output.txt" --taxid-paths "$DB/taxid-to-path.tsv" --output "$OUT/$RUN.calls.tsv"
    python -I "$TASKS/write_manifest.py" "$BENCH_ROOT" "$RUN" "$SETTING"
  done
done
