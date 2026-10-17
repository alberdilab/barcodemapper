"""Download one paired WGS ENA run and verify its published MD5 checksums."""

import argparse
import csv
import hashlib
import http.client
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen


def digest_file(path):
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
            md5.update(chunk)
            sha256.update(chunk)
            size += len(chunk)
    return size, md5.hexdigest(), sha256.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("run_accession")
    parser.add_argument("inputs_dir", type=Path)
    args = parser.parse_args()
    with args.report.open(newline="", encoding="utf-8") as handle:
        matches = [row for row in csv.DictReader(handle, delimiter="\t")
                   if row["run_accession"] == args.run_accession]
    if len(matches) != 1:
        parser.error("Run accession must occur exactly once in the ENA report")
    row = matches[0]
    if row["library_strategy"] != "WGS" or row["library_layout"] != "PAIRED":
        parser.error("Only paired WGS runs are accepted")
    urls = row["fastq_ftp"].split(";")
    md5s = row["fastq_md5"].split(";")
    sizes = row["fastq_bytes"].split(";")
    if not (len(urls) == len(md5s) == len(sizes) == 2):
        parser.error("ENA must provide exactly two FASTQ files, hashes and sizes")

    inputs_dir = args.inputs_dir.expanduser().resolve()
    raw_dir = inputs_dir / "raw"
    pair_dir = inputs_dir / "pairs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    pair_dir.mkdir(parents=True, exist_ok=True)
    manifest = inputs_dir / "raw" / f"{args.run_accession}.manifest.json"
    records = []
    for mate, (ftp, expected_md5, expected_size) in enumerate(zip(urls, md5s, sizes), 1):
        url = "https://" + ftp
        name = f"{args.run_accession}_{mate}.fastq.gz"
        if not ftp.endswith(name):
            parser.error(f"Unexpected ENA mate filename: {ftp}")
        target = raw_dir / name
        partial = raw_dir / (name + ".part")
        if target.exists():
            size, md5, sha256 = digest_file(target)
        else:
            for attempt in range(5):
                offset = partial.stat().st_size if partial.exists() else 0
                if offset == int(expected_size):
                    break
                if offset > int(expected_size):
                    raise ValueError(f"Partial FASTQ exceeds ENA size: {partial}")
                headers = {"User-Agent": "BarcodeMapper-benchmark/1.1.0"}
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                try:
                    with urlopen(Request(url, headers=headers), timeout=180) as response:
                        if offset and (response.status != 206 or
                                       not response.headers.get("Content-Range", "").startswith(f"bytes {offset}-")):
                            raise ValueError(f"ENA did not honor the resume range for {name}")
                        with partial.open("ab" if offset else "wb") as output:
                            for chunk in iter(lambda: response.read(8 * 1024 * 1024), b""):
                                output.write(chunk)
                    if partial.stat().st_size == int(expected_size):
                        break
                except (OSError, TimeoutError, http.client.IncompleteRead):
                    if attempt == 4:
                        raise
                time.sleep(min(30, 2 ** attempt))
            size, md5, sha256 = digest_file(partial)
        if size != int(expected_size) or md5 != expected_md5:
            raise ValueError(f"ENA checksum or size mismatch for {name}")
        if not target.exists():
            os.replace(str(partial), str(target))
        link = pair_dir / name
        if not link.exists():
            link.symlink_to(os.path.relpath(str(target), str(pair_dir)))
        elif link.resolve() != target:
            raise ValueError(f"Unexpected existing paired FASTQ link: {link}")
        records.append({"mate": mate, "url": url, "path": str(target),
                        "bytes": size, "ena_md5": md5, "sha256": sha256})
        print(f"{name}: {size} bytes, MD5 verified, SHA-256 {sha256}", flush=True)

    if manifest.exists():
        prior = json.loads(manifest.read_text(encoding="utf-8"))
        if prior.get("run_accession") != args.run_accession or prior.get("files") != records:
            raise ValueError(f"Existing manifest conflicts with verified files: {manifest}")
        print(f"Already verified: {manifest}", flush=True)
        return
    manifest.write_text(json.dumps({
        "run_accession": args.run_accession,
        "study_accession": row["study_accession"],
        "sample_accession": row["sample_accession"],
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "ena_report_sha256": digest_file(args.report)[2],
        "files": records,
    }, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
