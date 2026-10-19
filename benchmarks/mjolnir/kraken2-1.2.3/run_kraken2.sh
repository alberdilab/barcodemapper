#!/usr/bin/env bash
set -euo pipefail
source "$BENCH_ROOT/tasks/1.2.3-kraken2/common.sh"
mkdir -p "$OUT"
raw="$OUT/$RUN.output.txt"; report="$OUT/$RUN.report.txt"; timing="$OUT/$RUN.timing.json"
for f in "$raw" "$report" "$timing" "$OUT/$RUN.command.txt" "$OUT/$RUN.kraken2.log"; do
  [ ! -e "$f" ] || { echo "output exists, not overwriting: $f" >&2; exit 4; }
done
echo "host $(hostname) job ${SLURM_JOB_ID:-} setting $SETTING run $RUN db $DBNAME"
kraken2 --version | head -n 1
sha256sum "$CONFIG"
cmd=(kraken2 --db "$DB" --paired --confidence "$CONF" --threads 16
  --report "$report" --output "$raw" "$READS/${RUN}_1.fq.gz" "$READS/${RUN}_2.fq.gz")
printf '%q ' "${cmd[@]}" > "$OUT/$RUN.command.txt"; echo >> "$OUT/$RUN.command.txt"
cat "$OUT/$RUN.command.txt"
"$TASKS/time_classifier.sh" --label classifier-only --cpus 16 \
  --input "$READS/${RUN}_1.fq.gz" --input "$READS/${RUN}_2.fq.gz" --index "$DB/hash.k2d" --index "$DB/opts.k2d" --index "$DB/taxo.k2d" \
  --output "$timing" -- "${cmd[@]}" 2> "$OUT/$RUN.kraken2.log"
cat "$OUT/$RUN.kraken2.log"
echo "done: $(wc -l < "$raw") output lines"
