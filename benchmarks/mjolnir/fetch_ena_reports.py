"""Freeze ENA run metadata for the two public mock-community studies."""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


STUDIES = ("PRJEB45431", "PRJNA725368")
FIELDS = (
    "run_accession",
    "study_accession",
    "sample_accession",
    "sample_title",
    "experiment_title",
    "library_strategy",
    "library_layout",
    "instrument_model",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "ena-reports.json"
    paths = [output_dir / f"ena-{study}.tsv" for study in STUDIES]
    if manifest_path.exists() or any(path.exists() for path in paths):
        parser.error("ENA reports already exist; preserve the frozen snapshot")

    records = []
    for study, path in zip(STUDIES, paths):
        query = urlencode({
            "accession": study,
            "result": "read_run",
            "fields": ",".join(FIELDS),
            "format": "tsv",
        })
        url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
        with urlopen(Request(url, headers={"User-Agent": "BarcodeMapper-benchmark/1.1.0"}), timeout=90) as response:
            body = response.read()
        lines = body.decode("utf-8").splitlines()
        if not lines or lines[0].split("\t") != list(FIELDS) or len(lines) < 2:
            raise ValueError(f"Unexpected or empty ENA report for {study}")
        path.write_bytes(body)
        records.append({
            "study": study,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "url": url,
            "rows": len(lines) - 1,
            "sha256": hashlib.sha256(body).hexdigest(),
            "path": str(path),
        })
    manifest_path.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    for record in records:
        print(f"{record['study']}: {record['rows']} runs, SHA-256 {record['sha256']}")


if __name__ == "__main__":
    main()
