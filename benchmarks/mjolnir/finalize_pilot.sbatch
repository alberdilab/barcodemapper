#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${BENCH_ENV:?Export the project Conda prefix before submission}"
: "${RUN:?Export the ENA run accession before submission}"
: "${REFERENCE:?Export the reference basename before submission}"
: "${CONTROLLER_JOBID:?Export the controller job ID before submission}"
: "${CONTROLLER_TASK:?Export the controller task directory name before submission}"

"$BENCH_ENV/bin/python" "$BENCH_ROOT/plans/finalize_pilot.py" \
  "$BENCH_ROOT" "$RUN" "$REFERENCE" "$CONTROLLER_JOBID" "$CONTROLLER_TASK"
