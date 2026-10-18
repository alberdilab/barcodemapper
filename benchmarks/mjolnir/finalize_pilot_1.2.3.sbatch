#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail
: "${BENCH_ROOT:?}"
: "${BENCH_ENV:?}"
: "${RUN:?}"
: "${REFERENCE:?}"
: "${CONTROLLER_JOBID:?}"
: "${CONTROLLER_TASK:?}"
: "${MODE:?}"

"$BENCH_ENV/bin/python" "$BENCH_ROOT/tasks/1.2.3-pilots/finalize_pilot_1.2.3.py" \
  "$BENCH_ROOT" "$RUN" "$REFERENCE" "$CONTROLLER_JOBID" "$CONTROLLER_TASK" "$MODE"
