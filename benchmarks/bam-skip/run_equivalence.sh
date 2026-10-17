#!/usr/bin/env bash
# Check that streaming SAM extraction (sam_extract.py) reproduces the pipeline's own extract
# TSV for every finished hit-cap run (all samples, all-hits and every cap), and the pysam
# extractor for non-default thresholds. Reads the existing BAMs; no Bowtie2 rerun.
# Usage: run_equivalence.sh   (inside the benchmark env, node-local TMPDIR)
set -euo pipefail
ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
WORK_ROOT=$ROOT/projects/05-hit-cap/work
HERE=$(cd "$(dirname "$0")" && pwd)
PKG=${PKG:-$ROOT/projects/04-scaling/overlay-maxhits/barcodemapper/workflow/scripts}
OUT=$ROOT/projects/05-hit-cap/results/bam-skip/equivalence_all.tsv
TMP=${TMPDIR:-/tmp}/bamequiv.${SLURM_JOB_ID:-$$}
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
printf 'sample\tsetting\tthresholds\ttruth_lines\tstream_lines\tidentical\n' > "$OUT"

for bam in "$WORK_ROOT"/*/*/bowtie/*.bam; do
  setting=$(basename "$(dirname "$(dirname "$bam")")")
  sample=$(basename "$(dirname "$(dirname "$(dirname "$bam")")")")
  truth=$(ls "$(dirname "$(dirname "$bam")")"/extract/*.tsv)
  samtools view -h "$bam" | python "$HERE/sam_extract.py" --output "$TMP/stream.tsv"
  LC_ALL=C sort "$truth" > "$TMP/truth.sorted"
  LC_ALL=C sort "$TMP/stream.tsv" > "$TMP/stream.sorted"
  same=yes; cmp -s "$TMP/truth.sorted" "$TMP/stream.sorted" || same=NO
  printf '%s\t%s\tdefault(3,80)\t%s\t%s\t%s\n' "$sample" "$setting" "$(wc -l < "$truth")" "$(wc -l < "$TMP/stream.tsv")" "$same" >> "$OUT"
done

# Non-default thresholds against the pysam extractor itself.
for spec in "animal all" "fungal all" "fungal k5" "ERR6001167 k200"; do
  set -- $spec
  bam=$(ls "$WORK_ROOT/$1/$2"/bowtie/*.bam)
  for th in "1 50" "5 100" "0 30"; do
    set -- $spec
    read -r edit cov <<< "$th"
    python "$PKG/extract_taxonomy.py" --input "$bam" --output "$TMP/pysam.tsv" --max-edit-percent "$edit" --min-covered "$cov"
    samtools view -h "$bam" | python "$HERE/sam_extract.py" --output "$TMP/stream.tsv" --max-edit-percent "$edit" --min-covered "$cov"
    LC_ALL=C sort "$TMP/pysam.tsv" > "$TMP/pysam.sorted"
    LC_ALL=C sort "$TMP/stream.tsv" > "$TMP/stream.sorted"
    same=yes; cmp -s "$TMP/pysam.sorted" "$TMP/stream.sorted" || same=NO
    printf '%s\t%s\t(%s,%s)\t%s\t%s\t%s\n' "$1" "$2" "$edit" "$cov" "$(wc -l < "$TMP/pysam.tsv")" "$(wc -l < "$TMP/stream.tsv")" "$same" >> "$OUT"
  done
done
cat "$OUT"
