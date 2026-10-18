#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# Creates the BarcodeMapper 1.2.3 environment for Project 03 (software/env-1.2.3).
# The 1.1.2 environment (software/env) is not touched.

set -euo pipefail
BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
project_root="$BENCH_ROOT/projects/03-method-comparison"
software="$project_root/software"
project_env="$software/env-1.2.3"
wheel="barcodemapper-1.2.3-py3-none-any.whl"
wheel_sha=c468352b9e852b408100a21ba4837b5664f462cdd470fce6bb78e929ddbe8464

test -s "$software/environment-1.2.3.yml"
test ! -e "$project_env"
echo "$wheel_sha  $software/$wheel" | sha256sum -c -
printf '%s  %s\n' "$wheel_sha" "$wheel" > "$software/wheel.1.2.3.sha256"

export BENCH_ROOT PROJECT_ENV="$project_env" SOFTWARE="$software" WHEEL="$wheel"
bash -ic '
  set -eu
  useconda
  conda env create --prefix "$PROJECT_ENV" --file "$SOFTWARE/environment-1.2.3.yml"
  conda activate "$PROJECT_ENV"
  python -m pip install --no-deps "$SOFTWARE/$WHEEL"
  conda list --explicit > "$SOFTWARE/conda-linux-64.1.2.3.explicit.txt"
  conda list --json > "$SOFTWARE/conda-packages.1.2.3.json"
  echo "barcodemapper: $(barcodemapper --version 2>&1)"
  echo "bowtie2: $(bowtie2 --version 2>&1 | sed -n 1,1p)"
  echo "samtools: $(samtools --version 2>&1 | sed -n 1,1p)"
  echo "fastp: $(fastp --version 2>&1)"
  echo "kraken2: $(kraken2 --version 2>&1 | sed -n 1,1p)"
  echo "blastn: $(blastn -version 2>&1 | sed -n 1,1p)"
  echo "vsearch: $(vsearch --version 2>&1 | sed -n 1,1p)"
  echo "minimap2: $(minimap2 --version 2>&1)"
  ln -s "$BENCH_ROOT/software/sacct-wrapper.sh" "$PROJECT_ENV/bin/sacct"
  barcodemapper --version 2>&1 | grep -q "1\.2\.3"
'

(cd "$software" && sha256sum environment-1.2.3.yml wheel.1.2.3.sha256 "$wheel" \
  conda-linux-64.1.2.3.explicit.txt conda-packages.1.2.3.json > SHA256SUMS.1.2.3)
cat "$software/SHA256SUMS.1.2.3"
