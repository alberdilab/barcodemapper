#!/usr/bin/env bash
set -euo pipefail
source "$BENCH_ROOT/tasks/1.2.3-comparators/common.sh"
OUT="$RES/$TOOL"; mkdir -p "$OUT"
TC="$TASKS/time_classifier.sh"
echo "host $(hostname) job ${SLURM_JOB_ID:-} tool $TOOL run $RUN ref $REFNAME"
$TOOL --version > "$WORK/.ver.$$" 2>&1 || true; head -n 1 "$WORK/.ver.$$"; rm -f "$WORK/.ver.$$"
sha256sum "$CONFIG"
for m in 1 2; do
  q="$QDIR/${RUN}_$m.fa"
  raw="$OUT/$RUN.R$m.raw.${TOOL/minimap2/paf}"
  [ "$TOOL" = vsearch ] && raw="$OUT/$RUN.R$m.raw.tsv"
  timing="$OUT/$RUN.R$m.timing.json"
  [ ! -e "$raw" ] && [ ! -e "$timing" ] || { echo "output exists, not overwriting: $raw" >&2; exit 4; }
  if [ "$TOOL" = vsearch ]; then
    cmd=(vsearch --usearch_global "$q" --db "$REF" --userout "$raw"
      --userfields query+target+id+alnlen+mism+opens+gaps+ql+qlo+qhi+tl+caln+qcov
      --id 0.97 --iddef 2 --mincols 80 --strand both
      --maxaccepts 0 --maxrejects 1000 --top_hits_only
      --qmask none --dbmask none --fasta_width 0 --threads 16)
  else
    cmd=(minimap2 -x sr -t 16 -N 200 --secondary=yes -p 0.5 -I 32G -c -o "$raw" "$REF" "$q")
  fi
  printf '%q ' "${cmd[@]}" > "$OUT/$RUN.R$m.command.txt"; echo >> "$OUT/$RUN.R$m.command.txt"
  cat "$OUT/$RUN.R$m.command.txt"
  "$TC" --label classifier-only --cpus 16 --input "$READS/${RUN}_$m.fq.gz" --index "$REF" \
    --output "$timing" -- "${cmd[@]}"
  echo "R$m done: $(wc -l < "$raw") lines"
done
