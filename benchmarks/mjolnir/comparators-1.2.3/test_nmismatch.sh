#!/usr/bin/env bash
set -euo pipefail
P="$BENCH_ROOT/projects/03-method-comparison"; T="$BENCH_ROOT/tasks/1.2.3-comparators"
Wk="$P/work/comparators-1.2.3"; D="$Wk/test-nmismatch"; mkdir -p "$D"
REF="$BENCH_ROOT/references-1.2.3/default.fa.gz"
head -n 600 "$Wk/queries/ERR6001161_1.fa" > "$D/q300.fa"
vsearch --version 2>&1 | head -n 1 || true
vsearch --usearch_global "$D/q300.fa" --db "$REF" --userout "$D/raw.tsv" \
  --userfields query+target+id+alnlen+mism+opens+gaps+ql+qlo+qhi+tl+caln+qcov \
  --id 0.97 --iddef 2 --mincols 80 --strand both --maxaccepts 0 --maxrejects 1000 --top_hits_only \
  --n_mismatch --qmask none --dbmask none --fasta_width 0 --threads 16
python -I "$T/accept_hits.py" vsearch 1 "$D/raw.tsv" "$D/accepted.tsv"
zcat "$REF" | python -I "$T/vsearch_n_diagnostic.py" "$D/accepted.tsv" | tee "$D/diagnostic.txt"
# control: the same 300 queries in the first (invalid) run
python -I - "$P/results/pilots-1.2.3/vsearch/ERR6001161.R1.accepted.tsv" "$D/q300.fa" "$D/control_accepted.tsv" <<'PY'
import sys
ids = {l[1:].strip() for l in open(sys.argv[2]) if l[0] == ">"}
with open(sys.argv[3], "x") as o:
    for l in open(sys.argv[1]):
        if l.split("\t", 1)[0] in ids: o.write(l)
PY
zcat "$REF" | python -I "$T/vsearch_n_diagnostic.py" "$D/control_accepted.tsv" | tee "$D/control_diagnostic.txt"
