#!/usr/bin/env bash
# Run on the HPC: build a throwaway conda env from an env file, install this
# checkout into it, run the regression check, and print a one-block verdict.
#
# usage: mjolnir_dependency_check.sh REPO_DIR ENV_FILE [--keep]
# BM_SLURM_TEST=1 also repeats the check through Slurm (-s).
# BM_REGRESSION_ARGS adds options to regression_check.py (for example --update-golden).
set -uo pipefail

repo=${1:?usage: $0 REPO_DIR ENV_FILE [--keep]}
envfile=${2:?usage: $0 REPO_DIR ENV_FILE [--keep]}
keep=${3:-}
repo=$(cd "$repo" && pwd)
envfile=$(cd "$(dirname "$envfile")" && pwd)/$(basename "$envfile")

scratch=${BM_SCRATCH:-/projects/alberdilab/scratch/jpl786}
work=$(mktemp -d "$scratch/tmp_bm_depcheck.XXXXXX") || { echo "RESULT: FAIL (cannot create a work directory in $scratch)"; exit 4; }
env_dir=$work/env
status=1

cleanup() {
    if [ -n "$keep" ]; then
        echo "kept: $work"
    else
        rm -rf "$work"
    fi
}
trap cleanup EXIT

# conda is only set up in interactive shells on Mjolnir; load it directly.
# shellcheck disable=SC1091
source "${BM_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}" ||
    { echo "RESULT: FAIL (cannot load conda; set BM_CONDA_SH)"; exit 4; }

echo "== creating environment from $(basename "$envfile")"
if ! conda env create -q -p "$env_dir" -f "$envfile" >"$work/create.log" 2>&1; then
    tail -20 "$work/create.log"
    echo "RESULT: FAIL (environment could not be created; the versions may conflict)"
    exit 3
fi
if ! "$env_dir/bin/python" -m pip install -q --no-deps "$repo" >"$work/pip.log" 2>&1; then
    tail -20 "$work/pip.log"
    echo "RESULT: FAIL (barcodemapper did not install)"
    exit 3
fi

echo "== installed versions"
versions=$(conda list -p "$env_dir" 2>/dev/null |
    grep -E '^(python|snakemake|snakemake-executor-plugin-slurm|fastp|samtools|bowtie2|numpy|pandas|plotly|biopython|pysam|pyyaml) ' |
    awk '{printf "%s%s %s", sep, $1, $2; sep=", "}')
echo "$versions" | tr ',' '\n' | awk '{printf "  %-34s %s\n", $1, $2}'

run_check() {  # run_check LABEL [extra regression_check.py options]
    local label=$1; shift
    echo "== regression check ($label)"
    PATH="$env_dir/bin:$PATH" "$env_dir/bin/python" "$repo/scripts/regression_check.py" \
        --barcodemapper "$env_dir/bin/barcodemapper" "$@" ${BM_REGRESSION_ARGS:-}
}

local_result=PASS
run_check "local profile" || local_result=FAIL
slurm_result=skipped
if [ -n "${BM_SLURM_TEST:-}" ]; then
    slurm_result=PASS
    run_check "Slurm profile" --slurm || slurm_result=FAIL
fi

echo "RECORD|local=$local_result|slurm=$slurm_result|versions=$versions"
if [ "$local_result" = PASS ] && [ "$slurm_result" != FAIL ]; then
    echo "RESULT: PASS (results identical to the golden table; local $local_result, Slurm $slurm_result)"
    exit 0
fi
echo "RESULT: FAIL (local $local_result, Slurm $slurm_result; results differ or the pipeline failed, see above)"
exit 1
