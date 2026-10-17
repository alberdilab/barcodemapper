#!/usr/bin/env bash
# Benchmark the BAM step against streaming SAM straight into the extractor, reusing the
# existing hit-cap `-a` BAMs (no Bowtie2 rerun) and their extract TSVs as ground truth.
# Every variant reads the same SAM file, so only the post-Bowtie2 work differs:
#   bam_pysam : samtools view -b, then pysam extract_taxonomy.py   (current pipeline)
#   sam_py    : sam_extract.py on the SAM stream                   (proposed, no BAM)
#   sam_pysam : pysam reading SAM text                             (no BAM write)
#   sam_awk   : awk filter                                         (speed floor)
# Usage: run_bench.sh <animal|fungal>   env: REPS (3), THREADS (4)
# Run inside the benchmark env (samtools, pysam) with a node-local TMPDIR.
set -euo pipefail
SAMPLE=${1:?animal|fungal}
REPS=${REPS:-3}
THREADS=${THREADS:-4}
ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
W=$ROOT/projects/05-hit-cap/work/$SAMPLE/all
HERE=$(cd "$(dirname "$0")" && pwd)
export PKG=${PKG:-$ROOT/projects/04-scaling/overlay-maxhits/barcodemapper/workflow/scripts}
OUT=$ROOT/projects/05-hit-cap/results/bam-skip/$SAMPLE
WORK=${TMPDIR:-/tmp}/bamskip.${SLURM_JOB_ID:-$$}
mkdir -p "$OUT" "$WORK"
trap 'rm -rf "$WORK"' EXIT
BAM=$(ls "$W"/bowtie/*.bam | grep -v collated)
TRUTH=$(ls "$W"/extract/*.tsv)

samtools view -h -@ "$THREADS" "$BAM" > "$WORK/in.sam"
ALN=$(samtools view -c -@ "$THREADS" "$BAM")
printf 'alignments\t%s\nsam_bytes\t%s\nbam_bytes\t%s\n' "$ALN" "$(stat -c %s "$WORK/in.sam")" "$(stat -c %s "$BAM")" > "$OUT/sizes.tsv"
cat "$WORK/in.sam" > /dev/null   # warm the page cache

AWK='BEGIN{FS=OFS="\t"} /^@/{next} int($2/4)%2==1{next} $6=="*"{next}
{c=0; s=$6; while (match(s,/^[0-9]+[MIDNSHP=X]/)) {n=substr(s,1,RLENGTH-1)+0; o=substr(s,RLENGTH,1); if(o=="M"||o=="="||o=="X")c+=n; s=substr(s,RLENGTH+1)}
 if(c<80)next; nm=""; for(i=12;i<=NF;i++) if($i~/^NM:i:/){nm=substr($i,6)+0;break}
 if(nm==""||nm*100>3*c)next; split($3,p,"|"); print $1,p[1],p[2],p[3],nm,c}'

printf 'variant\trep\telapsed_s\tuser_s\tsys_s\tmax_rss_kb\n' > "$OUT/timings.tsv"
run() {  # variant rep command
  local v=$1 r=$2; shift 2
  python "$HERE/timeit.py" "$OUT/timings.tsv" "$v" "$r" "$*"
}
for rep in $(seq "$REPS"); do
  run bam_pysam "$rep" "samtools view -b -@ $THREADS $WORK/in.sam -o $WORK/x.bam && python $PKG/extract_taxonomy.py --input $WORK/x.bam --output $WORK/bam_pysam.tsv --max-edit-percent 3 --min-covered 80; rm -f $WORK/x.bam"
  run sam_py "$rep" "python $HERE/sam_extract.py --output $WORK/sam_py.tsv < $WORK/in.sam"
  run sam_pysam "$rep" "python $HERE/sam_pysam_extract.py --input $WORK/in.sam --output $WORK/sam_pysam.tsv"
  run sam_awk "$rep" "awk '$AWK' $WORK/in.sam > $WORK/sam_awk.tsv"
done

# Equivalence: every variant must equal the pipeline's own extract TSV as a multiset of lines.
LC_ALL=C sort "$TRUTH" > "$WORK/truth.sorted"
printf 'variant\tlines\tidentical_to_truth\n' > "$OUT/equivalence.tsv"
printf 'truth\t%s\tyes\n' "$(wc -l < "$TRUTH")" >> "$OUT/equivalence.tsv"
for v in bam_pysam sam_py sam_pysam sam_awk; do
  LC_ALL=C sort "$WORK/$v.tsv" > "$WORK/$v.sorted"
  if cmp -s "$WORK/truth.sorted" "$WORK/$v.sorted"; then same=yes; else same=NO; fi
  printf '%s\t%s\t%s\n' "$v" "$(wc -l < "$WORK/$v.tsv")" "$same" >> "$OUT/equivalence.tsv"
done
cat "$OUT/sizes.tsv" "$OUT/equivalence.tsv" "$OUT/timings.tsv"
