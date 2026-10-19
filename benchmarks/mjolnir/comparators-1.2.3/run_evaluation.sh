#!/usr/bin/env bash
# Stage 2: convert (committed converter) and evaluate the real-mock pilots.
# Run on the login node with BENCH_ROOT set; uses env-1.2.3 python. Refuses to overwrite.
set -euo pipefail
: "${BENCH_ROOT:?}"
P="$BENCH_ROOT/projects/03-method-comparison"; T="$BENCH_ROOT/tasks/1.2.3-comparators"
R="$P/results/pilots-1.2.3"; E="$R/evaluation"; W="$P/work/pilots-1.2.3"
PY="$P/software/env-1.2.3/bin/python"
CONV="$T/convert_calls.HEAD.py"; EVAL="$T/evaluate_calls.py"
EXPECTED="$BENCH_ROOT/projects/01-mock-communities/inputs/expected-taxa.tsv"
mkdir -p "$E/calls" "$E/hashes"
sha256sum "$CONV" "$EVAL" "$EXPECTED" > "$E/script_and_truth_hashes.txt"
declare -A REFNAME=( [SRR14351851]=all-fungi [ERR6001161]=default [ERR6001173]=default )
for RUN in SRR14351851 ERR6001161 ERR6001173; do
  REF="$BENCH_ROOT/references-1.2.3/${REFNAME[$RUN]}.fa.gz"
  FQ1="$W/$RUN-pilot-standard/fastp/${RUN}_1.fq.gz"; FQ2="$W/$RUN-pilot-standard/fastp/${RUN}_2.fq.gz"
  C="$E/calls"
  [ ! -e "$C/$RUN.barcodemapper_standard.tsv" ] || { echo "exists: evaluation already run for $RUN" >&2; exit 4; }
  "$PY" -I "$CONV" barcodemapper --read-ids "$FQ1" --input "$W/$RUN-pilot-standard/taxonomy/$RUN.tsv" --output "$C/$RUN.barcodemapper_standard.tsv"
  "$PY" -I "$CONV" barcodemapper --read-ids "$FQ1" --input "$W/$RUN-pilot-allhits/taxonomy/$RUN.tsv" --output "$C/$RUN.barcodemapper_allhits.tsv"
  "$PY" -I "$CONV" best-hit --read-ids "$FQ1" --input "$R/best-hit/$RUN.tsv" --output "$C/$RUN.best_hit.tsv"
  for pair in vsearch:vsearch:tsv:vsearch minimap2:minimap2:paf:minimap2_preset minimap2-sensitive:minimap2:paf:minimap2_sensitive; do
    IFS=: read dir sub ext name <<< "$pair"
    "$PY" -I "$CONV" $sub --read-ids "$FQ1" --input "$R/$dir/$RUN.accepted.$ext" \
      --accession-paths "$R/$dir/$RUN.hitmap.tsv" --output "$C/$RUN.$name.tsv"
    # converter-equivalence check against the earlier (working-tree converter) calls
    if cmp -s "$C/$RUN.$name.tsv" "$R/$dir/$RUN.calls.tsv"; then echo "$RUN $name: HEAD conversion identical to earlier calls.tsv"
    else echo "$RUN $name: HEAD conversion DIFFERS from earlier calls.tsv"; fi
  done
  # hashes manifest (real paths so evaluate_calls.py verifies the files)
  fq1s=$(sha256sum "$FQ1" | cut -d' ' -f1); fq2s=$(sha256sum "$FQ2" | cut -d' ' -f1); refs=$(sha256sum "$REF" | cut -d' ' -f1)
  H="$E/hashes/$RUN.hashes.tsv"
  printf 'method\tfastq_r1_sha256\tfastq_r2_sha256\treference_sha256\tfastq_r1_path\tfastq_r2_path\treference_path\n' > "$H"
  for name in barcodemapper_standard barcodemapper_allhits best_hit vsearch minimap2_preset minimap2_sensitive; do
    case $name in barcodemapper_allhits) d="$W/$RUN-pilot-allhits/fastp";; *) d="$W/$RUN-pilot-standard/fastp";; esac
    a=$(sha256sum "$d/${RUN}_1.fq.gz" | cut -d' ' -f1); b=$(sha256sum "$d/${RUN}_2.fq.gz" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $name $a $b $refs "$d/${RUN}_1.fq.gz" "$d/${RUN}_2.fq.gz" "$REF" >> "$H"
  done
  # truth for this run (blank has no expected taxa: a single placeholder row, so every call is unexpected)
  TRUTH="$EXPECTED"
  if [ "$RUN" = ERR6001173 ]; then
    TRUTH="$E/hashes/blank-expected-taxa.tsv"
    head -1 "$EXPECTED" > "$TRUTH"
    "$PY" -I - "$EXPECTED" >> "$TRUTH" <<'PY'
import csv, sys
h = next(csv.reader(open(sys.argv[1]), delimiter="\t"))
row = dict.fromkeys(h, "")
row.update(run_accession="ERR6001173", role="blank", expected_species="NO_EXPECTED_TAXA_BLANK_CONTROL",
           coverage_status="none", exact_species_reference_records="0", genus_reference_records="0", family_reference_records="0")
print("\t".join(row[k] for k in h))
PY
  fi
  args=(); for name in barcodemapper_standard barcodemapper_allhits best_hit vsearch minimap2_preset minimap2_sensitive; do args+=(--calls "$name=$C/$RUN.$name.tsv"); done
  for mf in 1 2; do
    "$PY" -I "$EVAL" "${args[@]}" --expected-taxa "$TRUTH" --run-accession "$RUN" --hashes-manifest "$H" --min-fragments $mf --out "$E/$RUN.minfrag$mf"
  done
done
"$PY" -I "$T/summarize_pilots.py" "$BENCH_ROOT"
