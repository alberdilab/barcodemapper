#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=controller-%j.out

set -euo pipefail
: "${BENCH_ENV:?Export an absolute BENCH_ENV Conda prefix before sbatch}"
if [[ "$BENCH_ENV" != /* ]]; then
    echo "BENCH_ENV must be an absolute path" >&2
    exit 2
fi
if (( $# == 0 )); then
    echo "Pass the command and its arguments after controller.sbatch" >&2
    exit 2
fi

# The Mjolnir useconda alias exists only in interactive Bash.
bash -ic 'useconda; conda activate "$BENCH_ENV"; exec "$@"' bash "$@"
