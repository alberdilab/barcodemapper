#!/usr/bin/env bash
# N-run diagnostic over all methods: one pass per reference. Read-only on inputs.
set -euo pipefail
: "${BENCH_ROOT:?}"
P="$BENCH_ROOT/projects/03-method-comparison"; R="$P/results/pilots-1.2.3"; N="$P/work/comparators-1.2.3/ndiag"
PY="$P/software/env-1.2.3/bin/python"; T="$BENCH_ROOT/tasks/1.2.3-comparators"
E="$R/evaluation/v2"; mkdir -p "$E"
spec() { run=$1; echo "$run.vsearch_invalid:vsearch:$R/vsearch/$run.accepted.tsv" "$run.vsearch_n_mismatch:vsearch:$R/vsearch-nmismatch/$run.accepted.tsv" \
  "$run.minimap2_preset:paf:$R/minimap2/$run.accepted.paf" "$run.minimap2_sensitive:paf:$R/minimap2-sensitive/$run.accepted.paf" \
  "$run.barcodemapper_standard:sam:$N/$run.bm.sam.txt"; }
zcat "$BENCH_ROOT/references-1.2.3/default.fa.gz" | "$PY" -I "$T/n_diag_all.py" "$E/n_run_diagnostic.default.tsv" $(spec ERR6001161) $(spec ERR6001173)
zcat "$BENCH_ROOT/references-1.2.3/all-fungi.fa.gz" | "$PY" -I "$T/n_diag_all.py" "$E/n_run_diagnostic.all-fungi.tsv" $(spec SRR14351851)
cat "$E/n_run_diagnostic.default.tsv" "$E/n_run_diagnostic.all-fungi.tsv" | awk 'NR==1||$1!="label"' > "$E/n_run_diagnostic.tsv"
cat "$E/n_run_diagnostic.tsv"
