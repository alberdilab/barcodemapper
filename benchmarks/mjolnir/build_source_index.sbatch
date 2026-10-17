#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"

python3 "$BENCH_ROOT/plans/build_source_index.py" \
  "$BENCH_ROOT/references/sources" \
  "$BENCH_ROOT/projects/02-reference-novelty/inputs"
