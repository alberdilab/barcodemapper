#!/usr/bin/env bash
# Print the Project 05 submission chain by default. CONFIRM_SUBMIT=yes runs it.
set -euo pipefail

BENCH_ROOT="${BENCH_ROOT:-/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks}"
PLAN_DIR="${PLAN_DIR:-$BENCH_ROOT/plans}"
MEM="${MEM:-1150G}"
CONFIRM_SUBMIT="${CONFIRM_SUBMIT:-}"
[[ "$CONFIRM_SUBMIT" == yes || -z "$CONFIRM_SUBMIT" ]] || {
    echo "CONFIRM_SUBMIT must be yes or empty" >&2
    exit 2
}
fungal_read_ids="${READ_IDS_SRR14351851:-}"
animal_read_ids="${READ_IDS_ERR6001161:-}"
blank_read_ids="${READ_IDS_ERR6001173:-}"
# Validate every export before the first sbatch call, so a bad subset cannot
# leave a partially submitted dependency chain.
for run in SRR14351851 ERR6001161 ERR6001173; do
    case "$run" in
        SRR14351851) read_ids="$fungal_read_ids" ;;
        ERR6001161) read_ids="$animal_read_ids" ;;
        ERR6001173) read_ids="$blank_read_ids" ;;
    esac
    if [[ "$read_ids" == *,* ]]; then
        echo "READ_IDS_$run contains a comma, which breaks sbatch --export" >&2
        exit 2
    fi
done
if [[ "$CONFIRM_SUBMIT" == yes ]]; then
    [[ -n "$fungal_read_ids" ]] || {
        echo "READ_IDS_SRR14351851 is required when CONFIRM_SUBMIT=yes" >&2; exit 2;
    }
    [[ -n "$animal_read_ids" ]] || {
        echo "READ_IDS_ERR6001161 is required when CONFIRM_SUBMIT=yes" >&2; exit 2;
    }
fi

submit() {
    if [[ "$CONFIRM_SUBMIT" == yes ]]; then
        command sbatch --parsable "$@"
    else
        printf 'sbatch --parsable'
        printf ' %q' "$@"
        printf '\n'
    fi
}

if [[ "$CONFIRM_SUBMIT" == yes ]]; then
    taxid_job=$(submit "$PLAN_DIR/05_taxid_paths.sbatch")
    taxid_job=${taxid_job%%;*}
    printf 'taxid_paths=%s\n' "$taxid_job"
else
    submit "$PLAN_DIR/05_taxid_paths.sbatch"
    taxid_job=PROJECT05_TAXID_JOB_ID
fi

previous_kraken=""
for run in SRR14351851 ERR6001161 ERR6001173; do
    case "$run" in
        SRR14351851) read_ids="$fungal_read_ids" ;;
        ERR6001161) read_ids="$animal_read_ids" ;;
        ERR6001173) read_ids="$blank_read_ids" ;;
    esac
    # Only the blank may use empty READ_IDS: search all retained pairs, but
    # cap it at the Project 03 subset size if large (not enforced here).
    submit --dependency="afterok:$taxid_job" \
        --export="ALL,RUN=$run,READ_IDS=$read_ids" \
        "$PLAN_DIR/05_blast_nt.sbatch"

    dependency="afterok:$taxid_job"
    if [[ -n "$previous_kraken" ]]; then
        dependency="${dependency}:$previous_kraken"
    fi
    if [[ "$CONFIRM_SUBMIT" == yes ]]; then
        previous_kraken=$(submit --dependency="$dependency" --mem="$MEM" \
            --export="ALL,RUN=$run,MEM=$MEM" \
            "$PLAN_DIR/05_kraken2_nt.sbatch")
        previous_kraken=${previous_kraken%%;*}
        printf 'kraken2_%s=%s\n' "$run" "$previous_kraken"
    else
        submit --dependency="$dependency" --mem="$MEM" \
            --export="ALL,RUN=$run,MEM=$MEM" \
            "$PLAN_DIR/05_kraken2_nt.sbatch"
        previous_kraken="PROJECT05_KRAKEN_${run}_JOB_ID"
    fi
done
