#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"

python3 "$BENCH_ROOT/plans/drop_empty_source_records.py" \
  "$BENCH_ROOT/references/sources" \
  "$BENCH_ROOT/references/cleaned-sources"
