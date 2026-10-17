#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2-01:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${CONTROLLER_JOBID:?Export the controller job ID before submission}"

project="$BENCH_ROOT/projects/04-scaling"
python3 "$BENCH_ROOT/plans/sample_scratch.py" \
  "$project/work/background-heavy/100k/rep1" \
  "$project/results/background-heavy/100k/rep1" \
  "$CONTROLLER_JOBID" --interval 30
