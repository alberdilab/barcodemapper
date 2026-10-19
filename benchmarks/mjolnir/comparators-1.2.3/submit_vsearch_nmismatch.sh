#!/usr/bin/env bash
set -euo pipefail
: "${BENCH_ROOT:?}"
T="$BENCH_ROOT/tasks/1.2.3-comparators"; L="$BENCH_ROOT/projects/03-method-comparison/logs"
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  r=$(sbatch --parsable --time=3-00:00:00 --job-name=cmp-vsnm-$RUN --output=$L/cmp-vsnm-$RUN-%j.out --export=ALL,RUN=$RUN,TOOL=vsearch $T/run_vsearch_nmismatch.sbatch)
  f=$(sbatch --parsable --dependency=afterok:$r --job-name=cmpfin-vsnm-$RUN --output=$L/cmpfin-vsnm-$RUN-%j.out --export=ALL,RUN=$RUN,TOOL=vsearch $T/finish_vsearch_nmismatch.sbatch)
  echo "$RUN search $r finish $f"
done
