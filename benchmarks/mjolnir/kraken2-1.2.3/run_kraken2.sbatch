#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=16
#SBATCH --mem=2G
#SBATCH --time=02:00:00
# Env: BENCH_ROOT, RUN, SETTING. DB is ~0.44 GB; 2G = DB x1.5 plus buffers.
set -euo pipefail
: "${BENCH_ROOT:?}" "${RUN:?}" "${SETTING:?}"
export PROJECT_ENV="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
exec bash -ic 'useconda; conda activate "$PROJECT_ENV"; exec bash "$BENCH_ROOT/tasks/1.2.3-kraken2/run_kraken2.sh"'
