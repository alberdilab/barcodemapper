#!/usr/bin/env bash
# Stage 2: convert (committed converter) and evaluate the real-mock pilots.
# Run on the login node with BENCH_ROOT set; uses env-1.2.3 python. Refuses to overwrite.
set -euo pipefail
: "${BENCH_ROOT:?}"
P="$BENCH_ROOT/projects/03-method-comparison"; T="$BENCH_ROOT/tasks/1.2.3-comparators"
R="$P/results/pilots-1.2.3"; E="$R/evaluation/v2"; C0="$R/evaluation/calls"; W="$P/work/pilots-1.2.3"
PY="$P/software/env-1.2.3/bin/python"
CONV="$T/convert_calls.HEAD.py"; EVAL="$T/evaluate_calls.py"
EXPECTED="$BENCH_ROOT/projects/01-mock-communities/inputs/expected-taxa.tsv"
mkdir -p "$E/hashes"
sha256sum "$CONV" "$EVAL" "$EXPECTED" > "$E/script_and_truth_hashes.v2.txt"
declare -A REFNAME=( [SRR14351851]=all-fungi [ERR6001161]=default [ERR6001173]=default )
# method key (calls file stem) -> label used in evaluator outputs
KEYS=(barcodemapper_standard barcodemapper_allhits best_hit vsearch_invalid vsearch_n_mismatch minimap2_preset minimap2_sensitive)
declare -A LABEL=( [barcodemapper_standard]=barcodemapper_standard [barcodemapper_allhits]=barcodemapper_allhits [best_hit]=best_hit
  [vsearch_invalid]="vsearch (invalid, N artifact)" [vsearch_n_mismatch]="vsearch (n_mismatch)"
  [minimap2_preset]=minimap2_preset [minimap2_sensitive]=minimap2_sensitive )
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  REF="$BENCH_ROOT/references-1.2.3/${REFNAME[$RUN]}.fa.gz"
  FQ1="$W/$RUN-pilot-standard/fastp/${RUN}_1.fq.gz"; FQ2="$W/$RUN-pilot-standard/fastp/${RUN}_2.fq.gz"
  [ ! -e "$E/$RUN.minfrag1" ] || { echo "exists: v2 evaluation already run for $RUN" >&2; exit 4; }
  # new conversion (HEAD converter) for the n_mismatch rerun only; all other calls are reused from evaluation/calls
  [ ! -e "$C0/$RUN.vsearch_n_mismatch.tsv" ] || { echo "exists: $C0/$RUN.vsearch_n_mismatch.tsv" >&2; exit 4; }
  "$PY" -I "$CONV" vsearch --read-ids "$FQ1" --input "$R/vsearch-nmismatch/$RUN.accepted.tsv" \
    --accession-paths "$R/vsearch-nmismatch/$RUN.hitmap.tsv" --output "$C0/$RUN.vsearch_n_mismatch.tsv"
  cmp -s "$C0/$RUN.vsearch_n_mismatch.tsv" "$R/vsearch-nmismatch/$RUN.calls.tsv" && echo "$RUN vsearch_n_mismatch: HEAD conversion identical to earlier calls.tsv" || echo "$RUN vsearch_n_mismatch: HEAD conversion DIFFERS from calls.tsv"
  cp -n "$C0/$RUN.vsearch_invalid.tsv" /dev/null 2>/dev/null || true
  [ -e "$C0/$RUN.vsearch_invalid.tsv" ] || cp "$C0/$RUN.vsearch.tsv" "$C0/$RUN.vsearch_invalid.tsv"
  refs=$(sha256sum "$REF" | cut -d' ' -f1)
  H="$E/hashes/$RUN.hashes.tsv"
  printf 'method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256\tfastq_r1_path\tfastq_r2_path\treference_path\n' > "$H"
  for key in "${KEYS[@]}"; do
    case $key in barcodemapper_allhits) d="$W/$RUN-pilot-allhits/fastp";; *) d="$W/$RUN-pilot-standard/fastp";; esac
    a=$(sha256sum "$d/${RUN}_1.fq.gz" | cut -d' ' -f1); b=$(sha256sum "$d/${RUN}_2.fq.gz" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${LABEL[$key]}" $a $b $refs "$d/${RUN}_1.fq.gz" "$d/${RUN}_2.fq.gz" "$REF" >> "$H"
  done
  TRUTH="$EXPECTED"
  if [ "$RUN" = ERR6001173 ]; then TRUTH="$R/evaluation/hashes/blank-expected-taxa.tsv"; fi   # placeholder row made by run_evaluation.sh
  args=(); for key in "${KEYS[@]}"; do args+=(--calls "${LABEL[$key]}=$C0/$RUN.$key.tsv"); done
  for mf in 1 2; do
    "$PY" -I "$EVAL" "${args[@]}" --expected-taxa "$TRUTH" --run-accession "$RUN" --hashes-manifest "$H" --min-fragments $mf --out "$E/$RUN.minfrag$mf"
  done
done
"$PY" -I "$T/summarize_pilots_v2.py" "$BENCH_ROOT"
