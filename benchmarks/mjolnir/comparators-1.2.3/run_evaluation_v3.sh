#!/usr/bin/env bash
# Final evaluation v3: v2 methods + Kraken2 (conf 0.0 and 0.1). Login node only; refuses to overwrite.
set -euo pipefail
: "${BENCH_ROOT:?}"
P="$BENCH_ROOT/projects/03-method-comparison"; T="$BENCH_ROOT/tasks/1.2.3-comparators"
R="$P/results/pilots-1.2.3"; E="$R/evaluation/v3"; C0="$R/evaluation/calls"; W="$P/work/pilots-1.2.3"
PY="$P/software/env-1.2.3/bin/python"; EVAL="$T/evaluate_calls.py"
EXPECTED="$BENCH_ROOT/projects/01-mock-communities/inputs/expected-taxa.tsv"
mkdir -p "$E/hashes"
sha256sum "$T/convert_calls.HEAD.py" "$EVAL" "$EXPECTED" "$T/run_evaluation_v3.sh" "$T/summarize_pilots_v3.py" > "$E/script_and_truth_hashes.v3.txt"
declare -A REFNAME=( [SRR14351851]=all-fungi [ERR6001161]=default [ERR6001173]=default )
KEYS=(barcodemapper_standard barcodemapper_allhits best_hit vsearch_n_mismatch minimap2_preset minimap2_sensitive kraken2_conf0.0 kraken2_conf0.1 vsearch_invalid)
declare -A LABEL=( [barcodemapper_standard]=barcodemapper_standard [barcodemapper_allhits]=barcodemapper_allhits [best_hit]=best_hit
  [vsearch_n_mismatch]="vsearch (n_mismatch)" [minimap2_preset]=minimap2_preset [minimap2_sensitive]=minimap2_sensitive
  [kraken2_conf0.0]="Kraken2 (conf 0.0)" [kraken2_conf0.1]="Kraken2 (conf 0.1)" [vsearch_invalid]="vsearch (invalid, N artifact)" )
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  REF="$BENCH_ROOT/references-1.2.3/${REFNAME[$RUN]}.fa.gz"
  [ ! -e "$E/$RUN.minfrag1" ] || { echo "exists: v3 already run for $RUN" >&2; exit 4; }
  declare -A CALLS=( [kraken2_conf0.0]="$R/kraken2/$RUN.calls.tsv" [kraken2_conf0.1]="$R/kraken2-conf0.1/$RUN.calls.tsv" )
  refs=$(sha256sum "$REF" | cut -d' ' -f1)
  H="$E/hashes/$RUN.hashes.tsv"
  printf 'method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256\tfastq_r1_path\tfastq_r2_path\treference_path\n' > "$H"
  args=()
  for key in "${KEYS[@]}"; do
    case $key in barcodemapper_allhits) d="$W/$RUN-pilot-allhits/fastp";; *) d="$W/$RUN-pilot-standard/fastp";; esac
    a=$(sha256sum "$d/${RUN}_1.fq.gz" | cut -d' ' -f1); b=$(sha256sum "$d/${RUN}_2.fq.gz" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${LABEL[$key]}" $a $b $refs "$d/${RUN}_1.fq.gz" "$d/${RUN}_2.fq.gz" "$REF" >> "$H"
    f="${CALLS[$key]:-$C0/$RUN.$key.tsv}"
    args+=(--calls "${LABEL[$key]}=$f")
  done
  TRUTH="$EXPECTED"; [ "$RUN" = ERR6001173 ] && TRUTH="$R/evaluation/hashes/blank-expected-taxa.tsv"
  for mf in 1 2; do
    "$PY" -I "$EVAL" "${args[@]}" --expected-taxa "$TRUTH" --run-accession "$RUN" --hashes-manifest "$H" --min-fragments $mf --out "$E/$RUN.minfrag$mf"
  done
  unset CALLS
done
"$PY" -I "$T/summarize_pilots_v3.py" "$BENCH_ROOT"
