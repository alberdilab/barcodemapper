"""1.2.3 variant (paths under pilots-1.2.3/). Assign each fastp-retained pair its best accepted BarcodeMapper alignment.

The prespecified proxy score is lowest NM, then longest covered bases. Exact
score ties resolve by accession, marker, and taxonomy in lexical order.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import sqlite3
from datetime import datetime, timezone
from itertools import groupby
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def accepted_groups(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for read_id, group in groupby(reader, key=lambda row: row[0]):
            rows = list(group)
            if any(len(row) != 6 for row in rows):
                raise ValueError(f"Invalid extracted alignment for {read_id}")
            candidates = []
            for row in rows:
                _, accession, marker, taxonomy, nm, covered = row
                if not accession or marker not in ("COI-5P", "ITS") or not taxonomy:
                    raise ValueError(f"Invalid reference metadata for {read_id}")
                candidates.append((int(nm), -int(covered), accession, marker, taxonomy))
            best = min(candidates)
            tied = [item for item in candidates if item[:2] == best[:2]]
            yield (read_id, best[4], best[2], best[3], best[0], -best[1],
                   len(tied), len({item[4] for item in tied}))


def fastq_ids(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError(f"Malformed FASTQ record in {path}")
            if len(sequence.strip()) != len(quality.strip()):
                raise ValueError(f"Sequence/quality length mismatch in {path}")
            read_id = header[1:].split()[0].removesuffix("/1")
            yield read_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bench_root", type=Path)
    parser.add_argument("run")
    args = parser.parse_args()
    root = args.bench_root.expanduser().resolve()
    work = root / "projects/03-method-comparison/work/pilots-1.2.3" / f"{args.run}-pilot-standard"
    extracted = work / "extract" / f"{args.run}.tsv"
    r1 = work / "fastp" / f"{args.run}_1.fq.gz"
    fastp_json = work / "fastp" / f"{args.run}.json"
    result_dir = root / "projects/03-method-comparison/results/pilots-1.2.3/best-hit"
    result_dir.mkdir(parents=True, exist_ok=True)
    output = result_dir / f"{args.run}.tsv"
    manifest = result_dir / f"{args.run}.manifest.json"
    database = result_dir / f"{args.run}.sqlite.partial"
    stage = result_dir / f"{args.run}.tsv.partial"
    if output.exists() or manifest.exists():
        parser.error("Ablation output already exists")
    if database.exists():
        database.unlink()
    if stage.exists():
        stage.unlink()
    connection = sqlite3.connect(str(database))
    hit_count = 0
    try:
        connection.execute("PRAGMA synchronous=NORMAL")
        connection.execute("CREATE TABLE hits (read_id TEXT PRIMARY KEY, taxonomy TEXT, accession TEXT, marker TEXT, nm INTEGER, covered INTEGER, score_ties INTEGER, taxonomy_ties INTEGER)")
        connection.execute("CREATE TABLE seen (read_id TEXT PRIMARY KEY)")
        for record in accepted_groups(extracted):
            connection.execute("INSERT INTO hits VALUES (?,?,?,?,?,?,?,?)", record)
            hit_count += 1
            if hit_count % 10000 == 0:
                connection.commit()
        connection.commit()
        pairs = assigned = ties = 0
        with stage.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("ReadID", "Taxonomy", "status", "Accession", "Marker", "NM",
                             "CoveredBases", "ScoreTieCount", "TaxonomyTieCount"))
            for read_id in fastq_ids(r1):
                connection.execute("INSERT INTO seen VALUES (?)", (read_id,))
                row = connection.execute(
                    "SELECT taxonomy,accession,marker,nm,covered,score_ties,taxonomy_ties FROM hits WHERE read_id=?",
                    (read_id,),
                ).fetchone()
                if row:
                    assigned += 1
                    ties += row[-1] > 1
                    writer.writerow((read_id, row[0], "accepted_alignment", *row[1:]))
                else:
                    writer.writerow((read_id, "", "no_accepted_alignment", "", "", "", "", "", ""))
                pairs += 1
        unobserved = connection.execute(
            "SELECT COUNT(*) FROM hits LEFT JOIN seen USING(read_id) WHERE seen.read_id IS NULL"
        ).fetchone()[0]
        if unobserved:
            raise ValueError(f"{unobserved} extracted read IDs not present in fastp R1")
        summary = json.loads(fastp_json.read_text(encoding="utf-8"))
        after = summary["summary"]["after_filtering"]["total_reads"]
        if after != pairs * 2:
            raise ValueError(f"fastp pair count mismatch: {pairs} pairs, {after} reads")
    finally:
        connection.close()
    database.unlink()
    os.replace(stage, output)
    manifest.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run": args.run,
        "score_policy": "lowest NM, then longest covered bases; lexical accession, marker, taxonomy for exact score ties",
        "input_pairs_after_fastp": pairs,
        "pairs_with_accepted_alignment": assigned,
        "pairs_without_accepted_alignment": pairs - assigned,
        "pairs_with_conflicting_taxonomy_score_tie": ties,
        "extracted_groups": hit_count,
        "sources": {str(path): sha256(path) for path in (extracted, r1, fastp_json)},
        "output_path": str(output),
        "output_sha256": sha256(output),
        "script_sha256": sha256(Path(__file__)),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"{args.run}: {assigned}/{pairs} pairs assigned by best accepted hit")


if __name__ == "__main__":
    main()
