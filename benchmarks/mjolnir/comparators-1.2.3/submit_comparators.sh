#!/usr/bin/env bash
# Submit prep, then per dataset and tool a search job and a finish job (afterok chain).
set -euo pipefail
: "${BENCH_ROOT:?}"
T="$BENCH_ROOT/tasks/1.2.3-comparators"
L="$BENCH_ROOT/projects/03-method-comparison/logs"
mkdir -p "$L"
prep=$(sbatch --parsable --output="$L/cmp-prep-%j.out" "$T/prep_queries.sbatch")
echo "prep $prep"
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  for TOOL in vsearch minimap2; do
    time=3-00:00:00; [ "$TOOL" = minimap2 ] && time=1-00:00:00
    run_id=$(sbatch --parsable --dependency=afterok:$prep --time=$time --job-name=cmp-$TOOL-$RUN \
      --output="$L/cmp-$TOOL-$RUN-%j.out" --export=ALL,RUN=$RUN,TOOL=$TOOL "$T/run_comparator.sbatch")
    fin_id=$(sbatch --parsable --dependency=afterok:$run_id --job-name=cmpfin-$TOOL-$RUN \
      --output="$L/cmpfin-$TOOL-$RUN-%j.out" --export=ALL,RUN=$RUN,TOOL=$TOOL "$T/finish_comparator.sbatch")
    echo "$RUN $TOOL search $run_id finish $fin_id"
  done
done
