"""Build a Kraken2 library and synthetic taxonomy from a BarcodeMapper FASTA."""

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path


RANKS = {rank: index for index, rank in enumerate("kpcofgs")}
RANK_NAMES = dict(zip(RANKS, ("kingdom", "phylum", "class", "order",
                              "family", "genus", "species")))
RANK_PATTERN = re.compile(r"^([kpcofgs])__(.*)$")
MARKERS = ("COI-5P", "ITS")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_records(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header = None
        sequence = []
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    if not sequence:
                        raise ValueError(f"Empty FASTA record before line {line_number}")
                    yield header, sequence
                header = line[1:]
                sequence = []
            elif line:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header at line {line_number}")
                sequence.append(line)
        if header is None:
            raise ValueError(f"Empty FASTA: {path}")
        if not sequence:
            raise ValueError(f"Empty last FASTA record: {path}")
        yield header, sequence


def reference_header(header):
    parts = header.split("|")
    if len(parts) != 3 or not parts[0] or any(char.isspace() for char in header):
        raise ValueError(f"Invalid BarcodeMapper reference header: {header}")
    accession, marker, taxonomy = parts
    if marker not in MARKERS:
        raise ValueError(f"Unsupported marker in reference header: {header}")
    previous_rank = -1
    for index, level in enumerate(taxonomy.split(";")):
        match = RANK_PATTERN.fullmatch(level)
        if match is None:
            raise ValueError(f"Invalid taxonomy level in reference header: {header}")
        rank, name = match.groups()
        rank_index = RANKS[rank]
        if (index == 0 and (rank != "k" or not name)) or rank_index <= previous_rank:
            raise ValueError(f"Invalid taxonomy rank order in reference header: {header}")
        previous_rank = rank_index
    return accession, marker, taxonomy


def source_inventory(fasta):
    counts = Counter()
    prefixes = set()
    for header, _ in fasta_records(fasta):
        _, marker, taxonomy = reference_header(header)
        counts[marker] += 1
        levels = taxonomy.split(";")
        prefixes.update(";".join(levels[:index]) for index in range(1, len(levels) + 1))
    return counts, {path: taxid for taxid, path in enumerate(sorted(prefixes), 2)}


def write_taxonomy(output, taxids):
    taxonomy_dir = output / "taxonomy"
    taxonomy_dir.mkdir(parents=True, exist_ok=True)
    with (taxonomy_dir / "nodes.dmp").open("w", encoding="utf-8") as nodes, \
            (taxonomy_dir / "names.dmp").open("w", encoding="utf-8") as names, \
            (output / "taxid-to-path.tsv").open("w", newline="", encoding="utf-8") as paths:
        writer = csv.writer(paths, delimiter="\t", lineterminator="\n")
        writer.writerow(("TaxID", "Taxonomy"))
        nodes.write("1\t|\t1\t|\tno rank\t|\n")
        names.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        writer.writerow((1, ""))
        for path, taxid in sorted(taxids.items(), key=lambda item: item[1]):
            parent = taxids.get(path.rpartition(";")[0], 1)
            level = path.rsplit(";", 1)[-1]
            rank = RANK_NAMES[level.split("__", 1)[0]]
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\n")
            names.write(f"{taxid}\t|\t{level}\t|\t\t|\tscientific name\t|\n")
            writer.writerow((taxid, path))


def write_library(fasta, output, taxids):
    with (output / "library.fa").open("w", encoding="utf-8") as library, \
            (output / "seqid2taxid.map").open("w", encoding="utf-8") as mapping:
        for index, (header, sequence) in enumerate(fasta_records(fasta), 1):
            accession, marker, taxonomy = reference_header(header)
            seqid = f"kraken:taxid|{taxids[taxonomy]}|{marker}-{index:09d}|{accession}"
            library.write(f">{seqid}\n")
            library.write("\n".join(sequence) + "\n")
            mapping.write(f"{seqid}\t{taxids[taxonomy]}\n")


def audit_outputs(source_counts, library, mapping):
    library_counts = Counter()
    map_counts = Counter()
    with mapping.open(encoding="utf-8") as handle:
        for record, row in zip_longest(fasta_records(library), handle):
            if record is None or row is None:
                raise ValueError("Kraken FASTA and seqid2taxid map record counts differ")
            header, _ = record
            parts = header.split("|")
            fields = row.rstrip("\n").split("\t")
            if (len(parts) != 4 or parts[0] != "kraken:taxid" or
                    len(fields) != 2 or fields != [header, parts[1]]):
                raise ValueError(f"Kraken FASTA and seqid2taxid map disagree: {header}")
            marker = parts[2].rsplit("-", 1)[0]
            if marker not in MARKERS:
                raise ValueError(f"Unexpected marker in Kraken sequence ID: {header}")
            library_counts[marker] += 1
            map_counts[marker] += 1
    if library_counts != source_counts or map_counts != source_counts:
        raise ValueError(f"Reference record count mismatch: source={dict(source_counts)}, "
                         f"library={dict(library_counts)}, map={dict(map_counts)}")
    return library_counts, map_counts


def build_inputs(fasta, output):
    fasta = Path(fasta).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    if not fasta.name.endswith(".gz"):
        raise ValueError("Expected a gzipped BarcodeMapper reference FASTA")
    source_counts, taxids = source_inventory(fasta)
    files = (output / "taxonomy/nodes.dmp", output / "taxonomy/names.dmp",
             output / "taxid-to-path.tsv", output / "library.fa",
             output / "seqid2taxid.map", output / "audit.json")
    for path in files:
        if path.exists():
            raise ValueError(f"Kraken input already exists: {path}")
    output.mkdir(parents=True, exist_ok=True)
    write_taxonomy(output, taxids)
    write_library(fasta, output, taxids)
    library_counts, map_counts = audit_outputs(source_counts, files[3], files[4])
    files[5].write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_fasta": str(fasta),
        "source_fasta_sha256": sha256(fasta),
        "source_records": sum(source_counts.values()),
        "source_records_by_marker": dict(sorted(source_counts.items())),
        "library_records": sum(library_counts.values()),
        "library_records_by_marker": dict(sorted(library_counts.items())),
        "map_records": sum(map_counts.values()),
        "map_records_by_marker": dict(sorted(map_counts.items())),
        "taxon_count": len(taxids) + 1,
        "outputs": {str(path.relative_to(output)): sha256(path) for path in files[:5]},
        "script_sha256": sha256(Path(__file__)),
    }, indent=2) + "\n", encoding="utf-8")
    return files[5]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(f"Built Kraken inputs: {build_inputs(args.fasta, args.output_dir)}")


if __name__ == "__main__":
    main()
