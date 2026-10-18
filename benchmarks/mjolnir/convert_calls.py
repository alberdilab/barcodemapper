"""Convert method calls to one ReadID, Taxonomy, status row per retained pair.

All commands take --read-ids (one ID per line, or an R1 .fq/.fastq[.gz]),
--input, and --output. Mapping TSVs have two columns, TaxID/Accession and
Taxonomy (a header is optional). BLASTn outfmt6 and vsearch userout must put
query and subject IDs in their first two columns; PAF uses its standard query
and target columns. Alignment inputs must already contain only accepted hits.
BLASTn uses bit scores and the tie/cap policy in blast_config.yaml.
"""

import argparse
import csv
import gzip
import math
import re
from pathlib import Path


def open_text(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8", newline="")


def pair_id(value):
    read_id = value.strip()
    if read_id.endswith(("/1", "/2")):
        read_id = read_id[:-2]
    if not read_id:
        raise ValueError("Empty read ID")
    return read_id


def retained_ids(path):
    fastq = path.suffix in (".fq", ".fastq") or path.suffix == ".gz" and path.with_suffix("").suffix in (".fq", ".fastq")
    with open_text(path) as handle:
        if fastq:
            while header := handle.readline():
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                if not header.startswith("@") or not plus.startswith("+") or not quality:
                    raise ValueError(f"Malformed FASTQ record in {path}")
                if len(sequence.strip()) != len(quality.strip()):
                    raise ValueError(f"Sequence/quality length mismatch in {path}")
                yield pair_id(header[1:].split()[0])
        else:
            for line_number, line in enumerate(handle, 1):
                if line_number == 1 and line.rstrip("\r\n") == "ReadID":
                    continue
                yield pair_id(line)


def resolved_prefix(taxonomy):
    if taxonomy.strip().lower() == "unclassified":
        return ""
    levels = []
    for level in taxonomy.strip().split(";"):
        if not level or level.endswith("__"):
            break
        levels.append(level)
    return ";".join(levels)


def common_taxonomy(first, second):
    shared = []
    for left, right in zip(first.split(";"), second.split(";")):
        if left != right:
            break
        shared.append(left)
    return ";".join(shared)


def add_call(calls, read_id, taxonomy, status):
    read_id = pair_id(read_id)
    taxonomy = resolved_prefix(taxonomy)
    previous = calls.get(read_id)
    if previous is None or previous[1] == "no_hit":
        calls[read_id] = (taxonomy, status)
    elif status != "no_hit":
        shared = common_taxonomy(previous[0], taxonomy) if previous[0] and taxonomy else ""
        calls[read_id] = (shared, "classified" if shared else "unclassified")


def load_paths(path, key_header):
    paths = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for number, row in enumerate(reader, 1):
            if number == 1 and len(row) == 2 and row[0].lower() == key_header.lower() and row[1].lower() in ("taxonomy", "path"):
                continue
            if (len(row) != 2 or not row[0] or
                    not resolved_prefix(row[1]) and not (key_header == "TaxID" and row == ["1", ""])):
                raise ValueError(f"Invalid mapping at {path}:{number}")
            if row[0] in paths and paths[row[0]] != resolved_prefix(row[1]):
                raise ValueError(f"Conflicting mapping for {row[0]!r}")
            paths[row[0]] = resolved_prefix(row[1])
    return paths


def table_calls(path, best_hit):
    calls = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ReadID", "Taxonomy"} | ({"status"} if best_hit else set())
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError(f"Missing columns in {path}: {sorted(required)}")
        for row in reader:
            if None in row.values():
                raise ValueError(f"Incomplete row in {path}")
            taxonomy = resolved_prefix(row["Taxonomy"])
            if best_hit:
                source_status = row["status"]
                if source_status not in ("accepted_alignment", "no_accepted_alignment"):
                    raise ValueError(f"Unknown best-hit status {source_status!r}")
                if source_status == "no_accepted_alignment" and taxonomy:
                    raise ValueError("No-hit row has taxonomy")
                status = "no_hit" if source_status == "no_accepted_alignment" else "classified" if taxonomy else "unclassified"
            else:
                status = "classified" if taxonomy else "unclassified"
            add_call(calls, row["ReadID"], taxonomy, status)
    return calls


def kraken_calls(path, paths):
    calls = {}
    with open_text(path) as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 3 or row[0] not in ("C", "U"):
                raise ValueError(f"Invalid Kraken2 row at {path}:{number}")
            taxid = row[2]
            if row[0] == "U" or taxid == "0":
                taxonomy, status = "", "unclassified"
            else:
                if taxid not in paths:
                    raise ValueError(f"Unmapped taxid {taxid!r} at {path}:{number}")
                taxonomy = paths[taxid]
                status = "classified" if taxonomy else "unclassified"
            add_call(calls, row[1], taxonomy, status)
    return calls


def alignment_calls(path, paths, format_name):
    calls = {}
    with open_text(path) as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < (12 if format_name == "minimap2" else 2):
                raise ValueError(f"Invalid {format_name} row at {path}:{number}")
            subject = row[5] if format_name == "minimap2" else row[1]
            accession = subject if subject in paths else subject.split("|", 1)[0]
            if accession not in paths:
                raise ValueError(f"Unmapped accession {subject!r} at {path}:{number}")
            add_call(calls, row[0], paths[accession], "classified")
    return calls


def blast_policy(path):
    settings = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\s*(max_target_seqs|score|ties|bit_score_tolerance):\s*([^#\s]+)", line)
        if match:
            settings[match[1]] = match[2]
    if settings.get("score") != "bit_score":
        raise ValueError(f"Unsupported BLAST score policy in {path}")
    ties = settings.get("ties")
    if ties == "exact_equal_score":
        tolerance = 0.0
    elif ties == "within_bit_score_tolerance":
        tolerance = float(settings.get("bit_score_tolerance", "nan"))
    else:
        raise ValueError(f"Unsupported BLAST tie policy in {path}")
    cap = int(settings.get("max_target_seqs", "0"))
    if not math.isfinite(tolerance) or tolerance < 0 or cap < 1:
        raise ValueError(f"Invalid BLAST tie tolerance or hit cap in {path}")
    return tolerance, cap


def blast_calls(path, paths, tolerance, cap):
    mate_hits = {}
    with open_text(path) as handle:
        for number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 12:
                raise ValueError(f"Invalid blastn row at {path}:{number}")
            subject = row[1]
            accession = subject if subject in paths else subject.split("|", 1)[0]
            if accession not in paths:
                raise ValueError(f"Unmapped accession {subject!r} at {path}:{number}")
            try:
                score = float(row[11])
            except ValueError as error:
                raise ValueError(f"Invalid BLAST bit score at {path}:{number}") from error
            if not math.isfinite(score):
                raise ValueError(f"Invalid BLAST bit score at {path}:{number}")
            read_id = pair_id(row[0])
            mate = row[0].strip()[-2:] if row[0].strip().endswith(("/1", "/2")) else ""
            mate_hits.setdefault((read_id, mate), []).append((score, subject, paths[accession]))
    selected = {}
    capped = set()
    for (read_id, _), hits in mate_hits.items():
        best = max(score for score, _, _ in hits)
        top = [(subject, taxonomy) for score, subject, taxonomy in hits
               if best - score <= tolerance]
        if len({subject for subject, _ in top}) >= cap:
            capped.add(read_id)
        else:
            selected.setdefault(read_id, []).extend(taxonomy for _, taxonomy in top)
    calls = {}
    for read_id, taxonomies in selected.items():
        if read_id not in capped:
            for taxonomy in taxonomies:
                add_call(calls, read_id, taxonomy, "classified")
    for read_id in capped:
        calls[read_id] = ("", "unclassified")
    return calls


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="method", required=True)
    for method in ("barcodemapper", "best-hit", "kraken2", "blastn", "vsearch", "minimap2"):
        command = commands.add_parser(method)
        command.add_argument("--read-ids", required=True, type=Path)
        command.add_argument("--input", required=True, type=Path)
        command.add_argument("--output", required=True, type=Path)
        if method == "kraken2":
            command.add_argument("--taxid-paths", required=True, type=Path)
        elif method in ("blastn", "vsearch", "minimap2"):
            command.add_argument("--accession-paths", required=True, type=Path)
            if method == "blastn":
                command.add_argument("--blast-config", type=Path,
                                     default=Path(__file__).with_name("blast_config.yaml"))
    args = parser.parse_args(argv)

    ids = list(retained_ids(args.read_ids))
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate retained read pair ID")
    if args.method in ("barcodemapper", "best-hit"):
        calls = table_calls(args.input, args.method == "best-hit")
    elif args.method == "kraken2":
        calls = kraken_calls(args.input, load_paths(args.taxid_paths, "TaxID"))
    elif args.method == "blastn":
        calls = blast_calls(args.input, load_paths(args.accession_paths, "Accession"),
                            *blast_policy(args.blast_config))
    else:
        calls = alignment_calls(args.input, load_paths(args.accession_paths, "Accession"), args.method)
    unexpected = calls.keys() - set(ids)
    if unexpected:
        raise ValueError(f"{len(unexpected)} call IDs not present in retained reads (e.g. {min(unexpected)!r})")

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("ReadID", "Taxonomy", "status"))
        for read_id in ids:
            taxonomy, status = calls.get(read_id, ("", "no_hit"))
            writer.writerow((read_id, taxonomy, status))


if __name__ == "__main__":
    main()
