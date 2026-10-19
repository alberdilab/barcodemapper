#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
# Env: BENCH_ROOT, RUN, TOOL. Exact acceptance filter, hit map, convert_calls.py,
# counts, hashes and manifest.
set -euo pipefail
: "${BENCH_ROOT:?}" "${RUN:?}" "${TOOL:?}"
export PROJECT_ENV="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
exec bash -ic 'useconda; conda activate "$PROJECT_ENV"; exec bash "$BENCH_ROOT/tasks/1.2.3-comparators/finish_minimap2_sensitive.sh"'
