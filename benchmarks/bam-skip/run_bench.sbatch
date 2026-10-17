#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --job-name=bm-bamskip
# Usage: sbatch --export=ALL,SAMPLE=animal run_bench.sbatch   (then fungal)
set -euo pipefail
ENV=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks/projects/04-scaling/software/env
export PATH="$ENV/bin:$PATH"
export TMPDIR=${TMPDIR:-/tmp}
bash "${SLURM_SUBMIT_DIR}/run_bench.sh" "${SAMPLE:?animal|fungal}"
