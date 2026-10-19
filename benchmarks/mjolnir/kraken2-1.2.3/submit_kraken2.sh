#!/usr/bin/env bash
# Submits the 6 classification jobs. Refuses to submit when outputs already exist.
set -euo pipefail
export BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
T="$BENCH_ROOT/tasks/1.2.3-kraken2"
for SETTING in kraken2 kraken2-conf0.1; do
  for RUN in SRR14351851 ERR6001161 ERR6001173; do
    [ ! -e "$BENCH_ROOT/projects/03-method-comparison/results/pilots-1.2.3/$SETTING/$RUN.output.txt" ] || { echo "exists: $SETTING $RUN" >&2; exit 4; }
    mkdir -p "$T/logs"
    sbatch --parsable --job-name="k2-$SETTING-$RUN" --output="$T/logs/%x-%j.out" \
      --export=ALL,RUN=$RUN,SETTING=$SETTING "$T/run_kraken2.sbatch"
  done
done
