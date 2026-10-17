#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail
: "${BENCH_ROOT:?Export the benchmark root before submission}"
: "${PROJECT:?Export the project name before submission}"

project_root="$BENCH_ROOT/projects/$PROJECT"
project_env="$project_root/software/env"
lock="$BENCH_ROOT/projects/01-mock-communities/software/conda-linux-64.explicit.txt"
mkdir -p "$project_root"/{software,inputs,work,results,logs}
cp "$BENCH_ROOT/plans/environment.yml" "$project_root/software/environment.yml"
cp "$BENCH_ROOT/plans/wheel.sha256" "$project_root/software/wheel.sha256"
cp "$BENCH_ROOT/software/barcodemapper-1.1.1-py3-none-any.whl" "$project_root/software/"
(cd "$project_root/software" && sha256sum -c wheel.sha256)
test -s "$lock"
test ! -e "$project_env"

export PROJECT_ROOT="$project_root" PROJECT_ENV="$project_env" BENCH_ROOT LOCK="$lock"
bash -ic '
  set -euo pipefail
  useconda
  conda create --yes --prefix "$PROJECT_ENV" --file "$LOCK"
  conda activate "$PROJECT_ENV"
  python -m pip install --no-deps "$PROJECT_ROOT/software/barcodemapper-1.1.1-py3-none-any.whl"
  conda list --explicit > "$PROJECT_ROOT/software/conda-linux-64.explicit.txt"
  conda list --json > "$PROJECT_ROOT/software/conda-packages.json"
  barcodemapper --version
'

(cd "$project_root/software" && sha256sum environment.yml conda-linux-64.explicit.txt conda-packages.json wheel.sha256 barcodemapper-1.1.1-py3-none-any.whl > SHA256SUMS)
