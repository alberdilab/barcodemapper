"""Make an auditable FASTA copy without zero-length source records.

BarcodeMapper v1.1.0 rejects zero-length records in its prepared reference.
The original downloaded FASTAs are retained unchanged.
"""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def transform(source, destination, excluded_handle):
    source_hash = hashlib.sha256()
    output_hash = hashlib.sha256()
    seen = kept = dropped = 0
    offset = 0
    record_offset = None
    header = None
    sequence_lines = []
    sequence_bases = 0

    def finish_record(output):
        nonlocal kept, dropped
        if header is None:
            return
        if sequence_bases:
            for item in (header, *sequence_lines):
                output.write(item)
                output_hash.update(item)
            kept += 1
        else:
            excluded_handle.write(f"{source.name}\t{record_offset}\t{header[1:].decode('utf-8').strip()}\n")
            dropped += 1

    with source.open("rb") as input_handle, destination.open("wb") as output_handle:
        for line in input_handle:
            source_hash.update(line)
            if line.startswith(b">"):
                finish_record(output_handle)
                header = line
                record_offset = offset
                sequence_lines = []
                sequence_bases = 0
                seen += 1
            elif header is None:
                if line.strip():
                    raise ValueError(f"Sequence before first header in {source}")
            else:
                sequence_lines.append(line)
                sequence_bases += len(line.strip())
            offset += len(line)
        finish_record(output_handle)
    return {
        "source_path": str(source),
        "source_bytes": source.stat().st_size,
        "source_sha256": source_hash.hexdigest(),
        "cleaned_bytes": destination.stat().st_size,
        "cleaned_sha256": output_hash.hexdigest(),
        "records_seen": seen,
        "records_kept": kept,
        "records_dropped_empty": dropped,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    sources = args.sources_dir.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    releases = json.loads((sources / "source-releases.json").read_text())
    summary_path = output / "empty-record-filter.json"
    excluded_path = output / "excluded-empty-records.tsv"
    if summary_path.exists() or excluded_path.exists():
        parser.error("Output manifest already exists")
    records = {}
    with (output / "excluded-empty-records.tsv.partial").open("w", encoding="utf-8") as excluded:
        excluded.write("source\tbyte_offset\theader\n")
        for name in ("bold", "unite"):
            source = sources / f"{name}.fasta"
            stage = output / f"{name}.fasta.partial"
            final = output / f"{name}.fasta"
            if final.exists():
                parser.error(f"Output already exists: {final}")
            result = transform(source, stage, excluded)
            expected = releases["sources"][name]["sha256"]
            if result["source_sha256"] != expected:
                raise ValueError(f"Frozen source checksum mismatch: {name}")
            os.replace(stage, final)
            result["cleaned_path"] = str(final)
            records[name] = result
            print(f"{name}: retained {result['records_kept']}, dropped {result['records_dropped_empty']} empty records", flush=True)
    os.replace(output / "excluded-empty-records.tsv.partial", excluded_path)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "drop_empty_source_records.py",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "operation": "Remove only FASTA records with no sequence bases; retain all bytes of every other record",
        "excluded_records_path": str(excluded_path),
        "excluded_records_sha256": hashlib.sha256(excluded_path.read_bytes()).hexdigest(),
        "sources": records,
    }
    summary_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    for path in (output / "bold.fasta", output / "unite.fasta", excluded_path, summary_path):
        path.chmod(0o444)


if __name__ == "__main__":
    main()
