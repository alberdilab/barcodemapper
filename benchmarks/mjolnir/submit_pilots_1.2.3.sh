#!/usr/bin/env bash
# BarcodeMapper 1.2.3 pilots: standard (-k 200 default) and --all-hits, three runs,
# each with a verify job; best-hit ablation after each primary verify.
# Never run twice (refuses if job-map.tsv exists). Set DRY_RUN=1 to echo only.
set -euo pipefail

BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
ENV123="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
PLAN="$BENCH_ROOT/plans"
TASK="$BENCH_ROOT/tasks/1.2.3-pilots"
REF="$BENCH_ROOT/references-1.2.3"
P01="$BENCH_ROOT/projects/01-mock-communities"
P03="$BENCH_ROOT/projects/03-method-comparison"
JOB_MAP="$TASK/job-map.tsv"
DRY_RUN=${DRY_RUN:-0}

FAKE=9000000
sbatch() {
    if [[ "$DRY_RUN" == 1 ]]; then FAKE=$((FAKE+1)); echo "sbatch $*" >&2; echo "$FAKE"
    else command sbatch "$@"; fi
}
if [[ "$DRY_RUN" == 1 ]]; then JOB_MAP=/dev/null
elif [[ -e "$JOB_MAP" ]]; then echo "already submitted: $JOB_MAP" >&2; exit 2
fi
test -x "$ENV123/bin/barcodemapper"; test -e "$ENV123/bin/sacct"
for r in SRR14351851 ERR6001161 ERR6001173; do
    test -s "$P01/inputs/pilot/${r}_1.fastq.gz"; test -s "$P01/inputs/pilot/${r}_2.fastq.gz"
done
if [[ "$DRY_RUN" != 1 ]]; then
    for l in fungal animal blank; do for m in standard allhits; do mkdir -p "$TASK/pilot-$l-$m"; done; done
    mkdir -p "$P03/results/pilots-1.2.3/best-hit" "$P03/work/pilots-1.2.3"
fi
printf 'stage\tjob_id\tdependency\tstatus\n' > "$JOB_MAP"
record() { printf '%s\t%s\t%s\tsubmitted\n' "$1" "$2" "$3" | tee -a "$JOB_MAP" >&2; }

submit_pilot() {   # run reference label mode
    local run="$1" name="$2" label="$3" mode="$4"
    local extra=() task="pilot-$label-$mode"
    [[ "$mode" == allhits ]] && extra=(--all-hits)
    local job
    job=$(sbatch --parsable \
        --export="ALL,BENCH_ENV=$ENV123" \
        --job-name="bm123-pilot-$label-$mode" \
        --output="$TASK/$task/controller-%j.out" \
        "$PLAN/controller.sbatch" \
        barcodemapper --slurm \
        -1 "$P01/inputs/pilot/${run}_1.fastq.gz" \
        -2 "$P01/inputs/pilot/${run}_2.fastq.gz" \
        -d "$REF/$name.fa.gz" \
        -o "$P03/results/pilots-1.2.3/$run-pilot-$mode.tsv" \
        -t "$P03/work/pilots-1.2.3/$run-pilot-$mode" \
        --standard --adapter-r1 auto --adapter-r2 auto ${extra[@]+"${extra[@]}"})
    record "pilot-$label-$mode" "$job" none
    local ver
    ver=$(sbatch --parsable --dependency="afterok:$job" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,BENCH_ENV=$ENV123,RUN=$run,REFERENCE=$name,CONTROLLER_JOBID=$job,CONTROLLER_TASK=$task,MODE=$mode" \
        --job-name="bm123-verify-$label-$mode" \
        --output="$TASK/verify-$label-$mode-%j.out" \
        "$TASK/finalize_pilot_1.2.3.sbatch")
    record "verify-$label-$mode" "$ver" "$job"
    LAST_VERIFY="$ver"
}

submit_besthit() {  # label run dependency
    local job
    job=$(sbatch --parsable --dependency="afterok:$3" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,BENCH_ENV=$ENV123,RUN=$2" \
        --job-name="bm123-besthit-$1" \
        --output="$TASK/besthit-$1-%j.out" \
        "$TASK/best_hit_ablation_1.2.3.sbatch")
    record "besthit-$1" "$job" "$3"
}

LAST_VERIFY=
for spec in "SRR14351851 all-fungi fungal" "ERR6001161 default animal" "ERR6001173 default blank"; do
    set -- $spec
    submit_pilot "$1" "$2" "$3" standard
    submit_besthit "$3" "$1" "$LAST_VERIFY"
    submit_pilot "$1" "$2" "$3" allhits
done
echo "job map: $JOB_MAP"
