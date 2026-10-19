#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=cmp-prep
# Verify fastp FASTQ hashes against the pilot manifests, check mate pairing,
# and write mate-labelled query FASTAs (">firstToken/1", ">firstToken/2").
set -euo pipefail
: "${BENCH_ROOT:?}"
mkdir -p "$BENCH_ROOT/projects/03-method-comparison/work/comparators-1.2.3/queries"
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  source "$BENCH_ROOT/tasks/1.2.3-comparators/common.sh"
  sums="$RES/$RUN-pilot-standard.SHA256SUMS"
  for m in 1 2; do
    fq="$READS/${RUN}_$m.fq.gz"
    got=$(sha256sum "$fq" | cut -d' ' -f1)
    want=$(awk -v f="$fq" '{p=$2; sub(/^\/maps/,"",p); g=f; if (p==g || $2==f) print $1}' "$sums" | head -1)
    [ -n "$want" ] || { echo "no manifest hash for $fq" >&2; exit 3; }
    [ "$got" = "$want" ] || { echo "HASH MISMATCH $fq $got $want" >&2; exit 3; }
    echo "$RUN R$m fastq sha256 $got == manifest"
    out="$QDIR/${RUN}_$m.fa"
    [ ! -e "$out" ] || { echo "exists, not overwriting: $out" >&2; exit 4; }
    zcat "$fq" | awk -v m="$m" 'NR%4==1{split($1,a," "); print ">" substr($1,2) "/" m} NR%4==2{print}' > "$out"
    echo "$RUN R$m fasta sha256 $(sha256sum "$out" | cut -d' ' -f1) records $(grep -c '^>' "$out")"
  done
  # mate pairing: identical ordered ID lists
  cmp <(grep '^>' "$QDIR/${RUN}_1.fa" | sed 's#/1$##') <(grep '^>' "$QDIR/${RUN}_2.fa" | sed 's#/2$##') \
    && echo "$RUN R1/R2 IDs identical and in the same order"
done
