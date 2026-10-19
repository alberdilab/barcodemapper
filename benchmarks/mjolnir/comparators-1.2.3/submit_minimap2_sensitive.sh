#!/usr/bin/env bash
set -euo pipefail
: "${BENCH_ROOT:?}"
T="$BENCH_ROOT/tasks/1.2.3-comparators"; L="$BENCH_ROOT/projects/03-method-comparison/logs"
export TOOL=minimap2
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  r=$(sbatch --parsable --time=1-00:00:00 --job-name=cmp-mm2sens-$RUN --output=$L/cmp-mm2sens-$RUN-%j.out --export=ALL,RUN=$RUN,TOOL=minimap2 $T/run_minimap2_sensitive.sbatch)
  f=$(sbatch --parsable --dependency=afterok:$r --job-name=cmpfin-mm2sens-$RUN --output=$L/cmpfin-mm2sens-$RUN-%j.out --export=ALL,RUN=$RUN,TOOL=minimap2 $T/finish_minimap2_sensitive.sbatch)
  echo "$RUN search $r finish $f"
done
