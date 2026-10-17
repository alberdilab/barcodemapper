#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --job-name=bm-bamequiv
# Usage: sbatch --chdir=<dir with run_equivalence.sh> run_equivalence.sbatch
set -euo pipefail
ENV=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks/projects/04-scaling/software/env
export PATH="$ENV/bin:$PATH"
export TMPDIR=${TMPDIR:-/tmp}
bash "${SLURM_SUBMIT_DIR}/run_equivalence.sh"
