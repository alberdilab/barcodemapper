#!/usr/bin/env bash
# Mjolnir starting points: check the compressed FASTA size and adjust memory
# (Bowtie2 indexing requests about 28 times that size) and wall time before use.
#SBATCH --job-name=bm-build-publish
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=3-00:00:00
#SBATCH --output=bm-build-publish-%j.out
#SBATCH --error=bm-build-publish-%j.err

set -euo pipefail
# Token values must stay out of shell traces, including when invoked with bash -x.
set +x

usage() {
    cat <<'EOF'
Usage: build_and_publish.sh --db-dir DIR [options]

Build the latest BOLD+UNITE database on Mjolnir and upload it to a Zenodo draft.
DIR must be under /projects/alberdilab/scratch/jpl786 for a live run
(override the scratch root with BM_SCRATCH_ROOT).

Options:
  --db-dir DIR       Parent of bmdb_<YYYYMMDD> folders (required)
  --env NAME_OR_PATH Conda environment (required for live runs unless BM_ENV is set);
                     shared environments: /projects/alberdilab/data/environments/conda/
  --record ID        Published Zenodo concept record ID for a new version
  --draft ID         Resume upload into an existing unpublished Zenodo draft
  --sandbox          Use the Zenodo sandbox
  --include-index    Upload the Bowtie2 index as well
  --publish          Publish immediately (default: leave a draft)
  --reuse DIR        Skip the build; upload this existing bmdb folder
  --dry-run          Build nothing and use the publisher's offline dry run;
                     needs --reuse or an existing matching folder in --db-dir
  --label NAME       Pass NAME to barcodemapper -n (builds only)
  -h, --help         Show this help

Set ZENODO_SANDBOX for --sandbox or ZENODO_PRODUCTION otherwise in
Codeome Settings > Secrets (enabled for HPC and this project). The script
reads the selected value from the environment, BM_SECRETS_FILE,
~/.codeome/secrets/barcodemapper-91de8232.env, or ~/.tokens, in that order.
Never pass a token as an argument. Set BM_CONDA_SH if Conda init is elsewhere.
EOF
}

fail() { printf 'Error: %s\n' "$*" >&2; exit 2; }

db_dir=
conda_env=${BM_ENV:-}
record=
draft=
reuse=
label=
sandbox=0
include_index=0
publish=0
dry_run=0

while (( $# )); do
    case $1 in
        --db-dir|--env|--record|--draft|--reuse|--label)
            option=$1
            (( $# >= 2 )) || fail "$option requires a value"
            [[ -n $2 && $2 != --* ]] || fail "$option requires a value"
            case $option in
                --db-dir) db_dir=$2 ;;
                --env) conda_env=$2 ;;
                --record) record=$2 ;;
                --draft) draft=$2 ;;
                --reuse) reuse=$2 ;;
                --label) label=$2 ;;
            esac
            shift 2 ;;
        --sandbox) sandbox=1; shift ;;
        --include-index) include_index=1; shift ;;
        --publish) publish=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) fail "unknown option: $1" ;;
    esac
done

[[ -n $db_dir ]] || { usage >&2; fail "--db-dir is required"; }
if (( ! dry_run )); then
    [[ -n $conda_env ]] || fail "--env is required unless BM_ENV is set"
fi
[[ -z $record || $record =~ ^[0-9]+$ ]] || fail "--record must contain digits only"
[[ -z $draft || $draft =~ ^[0-9]+$ ]] || fail "--draft must contain digits only"
[[ -z $record || -z $draft ]] || fail "--draft cannot be used with --record"
[[ -z $label || $label =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail "invalid --label"
[[ -z $reuse || -z $label ]] || fail "--label cannot be used with --reuse"
command -v python3 >/dev/null 2>&1 || fail "python3 is required"

# Resolve before checking the scratch boundary, including paths with '..'.
db_dir=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$db_dir")
if [[ -n $reuse ]]; then
    reuse=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$reuse")
    [[ -d $reuse && $(basename "$reuse") =~ ^bmdb_[0-9]{8}(_[A-Za-z0-9][A-Za-z0-9._-]*)?$ ]] ||
        fail "--reuse must point to an existing bmdb_<YYYYMMDD>[_label] folder"
fi

scratch=${BM_SCRATCH_ROOT:-/projects/alberdilab/scratch/jpl786}
# Resolve the root too: on Mjolnir /projects is a symlink to /maps/projects.
scratch=$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$scratch")
if (( ! dry_run )); then
    [[ $db_dir == "$scratch" || $db_dir == "$scratch/"* ]] ||
        fail "--db-dir must be under $scratch"
    if [[ -n $reuse ]]; then
        [[ $reuse == "$scratch/"* ]] || fail "--reuse must be under $scratch"
    fi
    if (( sandbox )); then
        token_name=ZENODO_SANDBOX
    else
        token_name=ZENODO_PRODUCTION
    fi
    unset token_value
    token_value=${!token_name:-}
    if [[ -z $token_value ]]; then
        # Read only the selected line; sourcing a file would export its other tokens.
        token_files=()
        [[ -z ${BM_SECRETS_FILE:-} ]] || token_files+=("$BM_SECRETS_FILE")
        token_files+=("$HOME/.codeome/secrets/barcodemapper-91de8232.env" "$HOME/.tokens")
        for token_file in "${token_files[@]}"; do
            [[ -f $token_file ]] || continue
            token_line=$(grep -E "^(export[[:space:]]+)?${token_name}=" "$token_file" | tail -n 1 || true)
            token_value=${token_line#*=}
            token_value=${token_value%\"}
            token_value=${token_value#\"}
            token_value=${token_value%\'}
            token_value=${token_value#\'}
            [[ -z $token_value ]] || break
        done
        unset token_files token_file token_line
    fi
    [[ -n $token_value ]] || fail "$token_name is required; set it in Codeome Settings > Secrets (enabled for HPC and this project)"
    unset ZENODO_SANDBOX ZENODO_PRODUCTION ZENODO_TOKEN
fi

if [[ -n ${BM_REPO:-} ]]; then
    repo=$BM_REPO
elif [[ -n ${SLURM_JOB_ID:-} ]]; then
    repo=${SLURM_SUBMIT_DIR:-}
    [[ -n $repo ]] || fail "SLURM_SUBMIT_DIR is required for a Slurm job; set BM_REPO"
else
    repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fi
publisher=$repo/scripts/publish_zenodo.py
[[ -f $publisher ]] || fail "publisher not found at $publisher; set BM_REPO to the repository root"

newest_database() {
    local folder name newest=
    for folder in "$db_dir"/bmdb_*; do
        [[ -d $folder ]] || continue
        name=${folder##*/}
        if [[ -n $label ]]; then
            [[ $name =~ ^bmdb_[0-9]{8}_(.+)$ && ${BASH_REMATCH[1]} == "$label" ]] || continue
        else
            [[ $name =~ ^bmdb_[0-9]{8}$ ]] || continue
        fi
        [[ -f $folder/$name.fa.gz && -f $folder/$name.manifest.json ]] || continue
        [[ -z $newest || $name > ${newest##*/} ]] && newest=$folder
    done
    [[ -n $newest ]] || return 1
    printf '%s\n' "$newest"
}

print_command() { printf '  '; printf '%q ' "$@"; printf '\n'; }

if (( ! dry_run )); then
    # The useconda alias is interactive-only on Mjolnir; load Conda directly.
    conda_sh=${BM_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}
    [[ -f $conda_sh ]] || fail "Conda init not found at $conda_sh; set BM_CONDA_SH"
    # shellcheck disable=SC1090
    set +u
    source "$conda_sh"
    command -v conda >/dev/null 2>&1 || fail "conda is unavailable after loading $conda_sh"
    conda activate "$conda_env"
    set -u
    command -v python >/dev/null 2>&1 || fail "python is missing from $conda_env"
    if [[ -z $reuse ]]; then
        for tool in barcodemapper snakemake bowtie2 bowtie2-build; do
            command -v "$tool" >/dev/null 2>&1 || fail "$tool is missing from $conda_env"
        done
        mkdir -p "$db_dir"
    fi
    python_cmd=python
else
    python_cmd=python3
fi

if [[ -n $reuse ]]; then
    database=$reuse
    printf 'Build: skipped (--reuse %s)\n' "$database"
else
    build_cmd=(barcodemapper --build --db-dir "$db_dir")
    [[ -z $label ]] || build_cmd+=(-n "$label")
    printf 'Build step:\n'
    print_command "${build_cmd[@]}"
    if (( dry_run )); then
        printf 'Build: skipped (--dry-run)\n'
    else
        # The CLI's default Snakemake work directory is relative to cwd.
        # Keep all large intermediate files under the same scratch directory.
        build_output=$(cd "$db_dir" && env -u ZENODO_SANDBOX -u ZENODO_PRODUCTION -u ZENODO_TOKEN "${build_cmd[@]}")
        printf '%s\n' "$build_output"
        built_database=
        while IFS= read -r line; do
            if [[ $line == 'Built latest database: '* ]]; then
                built_database=${line#'Built latest database: '}
            fi
        done <<< "$build_output"
    fi
    if [[ -n ${built_database:-} &&
          -f $built_database/$(basename "$built_database").manifest.json &&
          -f $built_database/$(basename "$built_database").fa.gz ]]; then
        database=$built_database
    else
        database=$(newest_database) || fail "no complete matching bmdb folder found in $db_dir; use --reuse for a specific folder"
    fi
fi

printf 'Database: %s\n' "$database"
publish_cmd=("$python_cmd" "$publisher" --database "$database")
[[ -z $record ]] || publish_cmd+=(--record "$record")
[[ -z $draft ]] || publish_cmd+=(--draft "$draft")
(( ! sandbox )) || publish_cmd+=(--sandbox)
(( ! include_index )) || publish_cmd+=(--include-index)
(( ! publish )) || publish_cmd+=(--publish)
(( ! dry_run )) || publish_cmd+=(--dry-run)
printf 'Publish step:\n'
print_command "${publish_cmd[@]}"
if (( dry_run )); then
    env -u ZENODO_SANDBOX -u ZENODO_PRODUCTION -u ZENODO_TOKEN "${publish_cmd[@]}"
else
    env -u ZENODO_SANDBOX -u ZENODO_PRODUCTION -u ZENODO_TOKEN ZENODO_TOKEN="$token_value" "${publish_cmd[@]}"
fi

if (( dry_run )); then
    printf '\nReview the planned files and source licences.\n'
    printf 'Run without --dry-run to create a draft; add --publish only after review.\n'
elif (( publish )); then
    printf '\nReview the published Zenodo record and source licences.\n'
else
    printf '\nReview the Zenodo draft and its source licences.\n'
    printf 'After review, publish that draft in the Zenodo UI, or rerun with --reuse %q --draft ID to resume its upload.\n' "$database"
fi
printf 'After publication, put the printed concept record ID in barcodemapper/zenodo.py DEFAULT_RECORD.\n'
