#!/usr/bin/env bash
# Rebuild the default and all-fungi references with BarcodeMapper 1.2.3 from the
# retained cleaned sources, into references-1.2.3/ (never references/).
set -euo pipefail

BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
ENV123="$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3"
PLAN="$BENCH_ROOT/plans"
TASK="$BENCH_ROOT/tasks/1.2.3-rebuild"
SRC="$BENCH_ROOT/references/cleaned-sources"
OUT="$BENCH_ROOT/references-1.2.3"
JOB_MAP="$TASK/job-map.tsv"

test -x "$ENV123/bin/barcodemapper"
test -e "$ENV123/bin/sacct"   # sacct wrapper link, required for snakemake status checks
mkdir -p "$OUT/audits"
# finalize_reference.py looks for cleaned-sources next to the FASTA to check provenance.
[[ -e "$OUT/cleaned-sources" ]] || ln -s "$SRC" "$OUT/cleaned-sources"

record() { printf '%s\t%s\t%s\tsubmitted\n' "$1" "$2" "$3" | tee -a "$JOB_MAP" >&2; }

for name in default all-fungi; do
    extra=()
    if [[ "$name" == all-fungi ]]; then extra=(-x k__Animalia -y k__Fungi); fi
    build=$(sbatch --parsable \
        --export="ALL,BENCH_ENV=$ENV123" \
        --job-name="bm123-build-$name" \
        --output="$TASK/build-$name-%j.out" \
        "$PLAN/controller.sbatch" \
        barcodemapper --build --slurm \
        -b "$SRC/bold.fasta" -u "$SRC/unite.fasta" \
        ${extra[@]+"${extra[@]}"} \
        -d "$OUT/$name.fa.gz" \
        -t "$OUT/work-$name")
    record "build-$name" "$build" none
    freeze=$(sbatch --parsable --dependency="afterok:$build" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,REFERENCE=$name" \
        --job-name="bm123-freeze-$name" \
        --output="$TASK/freeze-$name-%j.out" \
        "$TASK/finalize_reference_1.2.3.sbatch")
    record "freeze-$name" "$freeze" "$build"
    audit=$(sbatch --parsable --dependency="afterok:$freeze" \
        --export="ALL,BENCH_ROOT=$BENCH_ROOT,REFERENCE=$name" \
        --job-name="bm123-audit-$name" \
        --output="$TASK/audit-$name-%j.out" \
        "$TASK/audit_reference_1.2.3.sbatch")
    record "audit-$name" "$audit" "$freeze"
done
