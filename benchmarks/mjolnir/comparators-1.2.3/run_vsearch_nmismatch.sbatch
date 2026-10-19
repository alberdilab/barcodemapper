#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
# Env: BENCH_ROOT, RUN, TOOL (vsearch|minimap2). Searches R1 then R2 separately,
# each wrapped by time_classifier.sh. Never overwrites existing outputs.
set -euo pipefail
: "${BENCH_ROOT:?}" "${RUN:?}" "${TOOL:?}"
export PROJECT_ENV="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
exec bash -ic 'useconda; conda activate "$PROJECT_ENV"; exec bash "$BENCH_ROOT/tasks/1.2.3-comparators/run_vsearch_nmismatch.sh"'
