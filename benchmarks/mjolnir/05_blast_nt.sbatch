#!/usr/bin/env bash
#SBATCH --partition=cpuqueue
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=14-00:00:00
#SBATCH --output=05-blast-nt-%j.out
#SBATCH --error=05-blast-nt-%j.err
# The first run measures peak RSS; revise the 64G request if nt search needs more.

set -euo pipefail
: "${RUN:?Export the pilot run accession before submission}"

BENCH_ROOT="${BENCH_ROOT:-/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks}"
BENCH_ENV="${BENCH_ENV:-$BENCH_ROOT/projects/03-method-comparison/software/env-1.2.3}"
PLAN_DIR="${PLAN_DIR:-$BENCH_ROOT/plans}"
PROJECT_DIR="${PROJECT_DIR:-$BENCH_ROOT/projects/05-default-databases}"
OUT_DIR="${OUT_DIR:-$PROJECT_DIR/blast/$RUN}"
FASTP_ROOT="${FASTP_ROOT:-$BENCH_ROOT/projects/01-mock-communities/work}"
R1="${R1:-$FASTP_ROOT/$RUN-pilot-standard/fastp/${RUN}_1.fq.gz}"
R2="${R2:-$FASTP_ROOT/$RUN-pilot-standard/fastp/${RUN}_2.fq.gz}"
BLAST_DB="${BLAST_DB:-/maps/datasets/globe_databases/ncbi/nt/20250722/nt}"
BLAST_CONFIG="${BLAST_CONFIG:-$PLAN_DIR/blast_config.yaml}"
TAXID_PATHS="${TAXID_PATHS:-$PROJECT_DIR/taxonomy/blast-taxid-paths.tsv}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
PYTHON="$BENCH_ENV/bin/python"
BLASTN="$BENCH_ENV/bin/blastn"
BLAST_OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
export PATH="$BENCH_ENV/bin:$PATH"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "$THREADS" != "$SLURM_CPUS_PER_TASK" ]]; then
    echo "THREADS must equal SLURM_CPUS_PER_TASK for time_classifier.sh" >&2
    exit 2
fi
for input in "$PYTHON" "$BLASTN" "$PLAN_DIR/time_classifier.sh" \
             "$PLAN_DIR/build_nt_accession_paths.py" "$PLAN_DIR/convert_calls.py" \
             "$BLAST_CONFIG" "$TAXID_PATHS" "$R1" "$R2"; do
    [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }
done
[[ -f "${BLAST_DB}.000.nin" ]] || { echo "Missing BLAST nt index: ${BLAST_DB}.000.nin" >&2; exit 2; }
if [[ -n "${READ_IDS:-}" ]]; then
    [[ -f "$READ_IDS" ]] || { echo "Missing READ_IDS file: $READ_IDS" >&2; exit 2; }
fi
[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mkdir -p "$PROJECT_DIR" "$(dirname "$OUT_DIR")"
stage=$(mktemp -d "$(dirname "$OUT_DIR")/.blast-${RUN}.XXXXXXXX")
preserve_failed_stage() {
    local status=$? failed
    trap - EXIT
    set +e
    if [[ -d "$stage" ]]; then
        failed="${OUT_DIR}.failed.${SLURM_JOB_ID:-manual-$$}"
        if [[ ! -e "$failed" && ! -L "$failed" ]] && \
                mv -T -n "$stage" "$failed" && [[ ! -e "$stage" ]]; then
            echo "Job failed (exit $status); preserved outputs at $failed" >&2
        else
            echo "Job failed (exit $status); could not move stage, outputs remain at $stage" >&2
        fi
    fi
    exit "$status"
}
trap preserve_failed_stage EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

# Read the frozen, deliberately small YAML schema without adding a YAML package.
# Reject new keys/policies instead of silently running a different search.
"$PYTHON" - "$BLAST_CONFIG" "$stage/blast-settings.json" "$stage/blast-params.txt" \
    "$BLAST_DB" "$THREADS" "${READ_IDS:-}" "$BLAST_OUTFMT" <<'PY'
import hashlib
import json
import re
import sys
from pathlib import Path

source, output, params = map(Path, sys.argv[1:4])
database, threads, read_ids, outfmt = sys.argv[4:]
top = {}
policy = {}
in_policy = False
for line in source.read_text(encoding="utf-8").splitlines():
    line = line.split("#", 1)[0].rstrip()
    if not line.strip():
        continue
    if line == "lca_tie_policy:":
        in_policy = True
        continue
    match = re.fullmatch(r"(  )?([a-z_]+):\s*(\S+)", line)
    if not match:
        raise SystemExit(f"Unsupported BLAST config line: {line!r}")
    target = policy if match.group(1) else top
    if bool(match.group(1)) != in_policy:
        raise SystemExit(f"Unexpected BLAST config indentation: {line!r}")
    key, value = match.group(2), match.group(3).strip('"\'')
    if key in target:
        raise SystemExit(f"Duplicate BLAST config key: {key}")
    target[key] = value
if set(top) != {"blast_version", "task", "evalue", "min_percent_identity",
                "min_query_coverage_percent", "max_target_seqs"}:
    raise SystemExit(f"Unexpected BLAST config keys: {sorted(top)}")
expected_policy = {
    "score": "bit_score", "ties": "exact_equal_score",
    "paired_reads": "union_top_scoring_accepted_hits_per_mate",
    "assignment": "lowest_common_ancestor_of_tied_taxonomy_paths",
    "accepted_assignment": "classified", "no_accepted_hits": "no_hit",
    "top_tie_reaches_max_target_seqs": "unclassified",
}
if policy != expected_policy or top["task"] != "blastn" or top["blast_version"] != "2.16.0":
    raise SystemExit("Unsupported BLAST version, task, or LCA policy")
try:
    evalue = float(top["evalue"])
    identity = float(top["min_percent_identity"])
    coverage = float(top["min_query_coverage_percent"])
    cap = int(top["max_target_seqs"])
except ValueError as error:
    raise SystemExit(f"Invalid BLAST threshold: {error}") from error
if not (0 < evalue and 0 <= identity <= 100 and 0 <= coverage <= 100 and cap > 0):
    raise SystemExit("BLAST thresholds out of range")
values = [top["evalue"], top["min_percent_identity"],
          top["min_query_coverage_percent"], top["max_target_seqs"]]
params.write_text("\n".join(values) + "\n", encoding="utf-8")
output.write_text(json.dumps({
    "blast_config_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "blast_version": top["blast_version"], "task": top["task"],
    "evalue": evalue, "perc_identity": identity,
    "qcov_hsp_perc": coverage, "max_target_seqs": cap,
    "database": database, "num_threads": int(threads),
    "read_ids_file": read_ids or None,
    "outfmt": outfmt,
    "lca_tie_policy": policy,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY
mapfile -t blast_params < "$stage/blast-params.txt"
rm "$stage/blast-params.txt"
version=$("$BLASTN" -version)
[[ "$version" == *"2.16.0"* ]] || { echo "Expected blastn 2.16.0, got $version" >&2; exit 2; }

# This writes one FASTA record per mate and one ReadID per retained pair.
# Selection is by pair ID, never by an independently filtered mate stream.
"$PYTHON" - "$R1" "$R2" "${READ_IDS:-}" "$stage/query.fa" "$stage/read-ids.tsv" <<'PY'
import gzip
import sys
from itertools import zip_longest
from pathlib import Path

r1, r2, ids_file, fasta_file, retained_file = sys.argv[1:]
requested = None
if ids_file:
    requested = set()
    with Path(ids_file).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            value = line.strip()
            if number == 1 and value == "ReadID":
                continue
            if value.endswith(("/1", "/2")):
                value = value[:-2]
            if not value or value in requested:
                raise SystemExit(f"Empty or duplicate READ_IDS entry at line {number}")
            requested.add(value)

def records(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while header := handle.readline():
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            if not header.startswith("@") or not plus.startswith("+") or not qual:
                raise SystemExit(f"Malformed FASTQ record in {path}")
            seq, qual = seq.rstrip("\r\n"), qual.rstrip("\r\n")
            if len(seq) != len(qual):
                raise SystemExit(f"Sequence/quality length mismatch in {path}")
            identifier = header[1:].split()[0]
            if identifier.endswith(("/1", "/2")):
                identifier = identifier[:-2]
            yield identifier, seq

seen = set()
with Path(fasta_file).open("w", encoding="utf-8") as fasta, \
        Path(retained_file).open("w", encoding="utf-8") as retained:
    for pair in zip_longest(records(r1), records(r2)):
        if pair[0] is None or pair[1] is None or pair[0][0] != pair[1][0]:
            raise SystemExit("FASTQ mates are missing or out of order")
        identifier = pair[0][0]
        if identifier in seen:
            raise SystemExit(f"Duplicate FASTQ pair ID: {identifier}")
        seen.add(identifier)
        if requested is not None and identifier not in requested:
            continue
        retained.write(identifier + "\n")
        fasta.write(f">{identifier}/1\n{pair[0][1]}\n>{identifier}/2\n{pair[1][1]}\n")
if requested is not None and requested - seen:
    raise SystemExit(f"{len(requested - seen)} READ_IDS entries are absent from FASTQs")
if not seen or not Path(retained_file).stat().st_size:
    raise SystemExit("No read pairs selected")
PY

sha256sum "$R1" "$R2" "$BLAST_CONFIG" "$TAXID_PATHS" > "$stage/input-sha256.txt"
if [[ -n "${READ_IDS:-}" ]]; then sha256sum "$READ_IDS" >> "$stage/input-sha256.txt"; fi
"$PLAN_DIR/time_classifier.sh" --label classifier-only --cpus "$THREADS" \
    --input "$stage/query.fa" --index "$(dirname "$BLAST_DB")" \
    --output "$stage/classifier-timing.json" -- \
    "$BLASTN" -task blastn -db "$BLAST_DB" -query "$stage/query.fa" \
    -evalue "${blast_params[0]}" -perc_identity "${blast_params[1]}" \
    -qcov_hsp_perc "${blast_params[2]}" -max_target_seqs "${blast_params[3]}" \
    -num_threads "$THREADS" \
    -outfmt "$BLAST_OUTFMT" \
    -out "$stage/blast-hits.tsv"
"$PYTHON" "$PLAN_DIR/build_nt_accession_paths.py" \
    --blast "$stage/blast-hits.tsv" --taxid-paths "$TAXID_PATHS" \
    --output "$stage/accession-paths.tsv" --audit "$stage/accession-paths.audit.json" \
    --unmapped-as-empty
"$PYTHON" "$PLAN_DIR/convert_calls.py" blastn \
    --read-ids "$stage/read-ids.tsv" --input "$stage/blast-hits.tsv" \
    --accession-paths "$stage/accession-paths.tsv" --blast-config "$BLAST_CONFIG" \
    --output "$stage/calls.tsv"

[[ ! -e "$OUT_DIR" && ! -L "$OUT_DIR" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
mv -T -n "$stage" "$OUT_DIR"
[[ ! -e "$stage" ]] || { echo "Refusing to overwrite $OUT_DIR" >&2; exit 2; }
trap - EXIT
