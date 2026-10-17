"""Inventory frozen marker references for mock coverage and holdout design."""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


RANKS = ("k", "p", "c", "o", "f", "g", "s")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    fasta = args.fasta.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths_file = output_dir / "taxonomy-paths.tsv"
    species_file = output_dir / "species.tsv"
    summary_file = output_dir / "summary.json"
    if any(path.exists() for path in (paths_file, species_file, summary_file)):
        parser.error("Reference audit already exists")
    base = fasta.with_suffix("").with_suffix("")
    snapshot = json.loads(Path(str(base) + ".snapshot.json").read_text(encoding="utf-8"))
    fasta_sha = next(item["sha256"] for item in snapshot["artifacts"]
                     if item["path"] == str(fasta))
    counts = Counter()
    species = Counter()
    marker_counts = Counter()
    blank_ranks = Counter()
    with gzip.open(str(fasta), "rb") as handle:
        for line in handle:
            if not line.startswith(b">"):
                continue
            header = line[1:].decode("utf-8").strip()
            accession, marker, taxonomy = header.split("|", 2)
            if not accession or marker not in ("COI-5P", "ITS"):
                raise ValueError(f"Unexpected reference header: {header}")
            levels = taxonomy.split(";")
            rank_map = {}
            for level in levels:
                rank, name = level.split("__", 1)
                if rank not in RANKS:
                    raise ValueError(f"Invalid rank: {level}")
                rank_map[rank] = name
                if not name:
                    blank_ranks[(marker, rank)] += 1
            counts[(marker, taxonomy)] += 1
            marker_counts[marker] += 1
            if rank_map.get("s"):
                species[(marker, rank_map["s"], rank_map.get("g", ""))] += 1
    with paths_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Marker", "Taxonomy", "Records"))
        for (marker, taxonomy), count in sorted(counts.items()):
            writer.writerow((marker, taxonomy, count))
    with species_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("Marker", "Species", "Genus", "Records"))
        for (marker, name, genus), count in sorted(species.items()):
            writer.writerow((marker, name, genus, count))
    summary_file.write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "fasta": str(fasta),
        "fasta_sha256": fasta_sha,
        "records_by_marker": dict(sorted(marker_counts.items())),
        "unique_taxonomy_paths": len(counts),
        "species_labels": len(species),
        "blank_ranks": {f"{marker}:{rank}": value
                        for (marker, rank), value in sorted(blank_ranks.items())},
        "taxonomy_paths_sha256": sha256(paths_file),
        "species_sha256": sha256(species_file),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Inventoried {sum(marker_counts.values())} records from {fasta}")


if __name__ == "__main__":
    main()
