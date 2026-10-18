"""Evaluate converted pair calls against synthetic truth or real mock members.

Example: python evaluate_calls.py --calls barcodemapper=calls.tsv --truth truth.tsv --out results

Repeat --calls for each method. Synthetic truth is a TSV with ReadID, Taxonomy,
and off_target (yes/no); off-target rows have empty Taxonomy. Mock evaluation
uses the expected-taxa.tsv made by annotate_expected_taxa.py and requires
--run-accession. The optional hashes manifest is a TSV with method,
fastq_r1_sha256, fastq_r2_sha256, and reference_sha256, one row per method.
Optional matching *_path columns make the listed file hashes mandatory to verify.
All inputs must describe the same retained read pairs.
"""

import argparse
import csv
import hashlib
import re
from collections import Counter
from pathlib import Path


RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
PREFIXES = dict(zip("kpcofgs", RANKS))
EVALUATED_RANKS = RANKS[1:]
MOCK_RANKS = ("family", "genus", "species")
CALL_FIELDS = ("ReadID", "Taxonomy", "status")
HASH_FIELDS = ("fastq_r1_sha256", "fastq_r2_sha256", "reference_sha256")
HASH_PATH_FIELDS = ("fastq_r1_path", "fastq_r2_path", "reference_path")


def read_tsv(path, required):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"Missing columns in {path}: {sorted(set(required) - set(reader.fieldnames or ()))}")
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"Incomplete row in {path}")
    return rows


def taxonomy_levels(value):
    if not value or value == "Unclassified":
        return {}
    levels = {}
    last_index = -1
    unresolved = False
    for part in value.split(";"):
        if "__" not in part:
            raise ValueError(f"Invalid taxonomy level {part!r}")
        prefix, name = part.split("__", 1)
        if prefix not in PREFIXES:
            raise ValueError(f"Invalid taxonomy rank {prefix!r}")
        rank = PREFIXES[prefix]
        index = RANKS.index(rank)
        if index <= last_index or name and unresolved:
            raise ValueError(f"Out-of-order taxonomy {value!r}")
        last_index = index
        if name:
            levels[rank] = name
        else:
            unresolved = True
    return levels


def taxon_key(levels, rank):
    return tuple((item, levels[item]) for item in RANKS[:RANKS.index(rank) + 1]
                 if item in levels)


def load_calls(path):
    calls = {}
    for row in read_tsv(path, CALL_FIELDS):
        read_id = row["ReadID"]
        if not read_id or read_id in calls:
            raise ValueError(f"Empty or duplicate ReadID {read_id!r} in {path}")
        status = row["status"]
        if status not in ("classified", "unclassified", "no_hit"):
            raise ValueError(f"Invalid status {status!r} in {path}")
        levels = taxonomy_levels(row["Taxonomy"])
        if bool(levels) != (status == "classified"):
            raise ValueError(f"Taxonomy/status mismatch for {read_id!r} in {path}")
        calls[read_id] = (levels, status)
    return calls


def load_truth(path):
    truth = {}
    for row in read_tsv(path, ("ReadID", "Taxonomy", "off_target")):
        read_id = row["ReadID"]
        if not read_id or read_id in truth:
            raise ValueError(f"Empty or duplicate ReadID {read_id!r} in {path}")
        flag = row["off_target"].strip().lower()
        if flag not in ("yes", "no", "true", "false", "1", "0"):
            raise ValueError(f"Invalid off_target flag {flag!r} for {read_id!r}")
        off_target = flag in ("yes", "true", "1")
        levels = taxonomy_levels(row["Taxonomy"])
        if bool(levels) == off_target:
            raise ValueError(f"Truth taxonomy/off_target mismatch for {read_id!r}")
        truth[read_id] = (levels, off_target)
    if not truth:
        raise ValueError("Truth has no read pairs")
    return truth


def check_hashes(path, methods):
    rows = read_tsv(path, ("method", *HASH_FIELDS))
    hashes = {}
    actual_hashes = {}
    for row in rows:
        method = row["method"]
        if not method or method in hashes:
            raise ValueError(f"Empty or duplicate method {method!r} in {path}")
        values = tuple(row[field].lower() for field in HASH_FIELDS)
        if any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in values):
            raise ValueError(f"Invalid sha256 for {method!r} in {path}")
        for field, path_field, expected in zip(HASH_FIELDS, HASH_PATH_FIELDS, values):
            if path_field not in row or not row[path_field]:
                continue
            file_path = Path(row[path_field])
            if not file_path.is_absolute():
                file_path = path.parent / file_path
            if file_path not in actual_hashes:
                digest = hashlib.sha256()
                with file_path.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                        digest.update(chunk)
                actual_hashes[file_path] = digest.hexdigest()
            if actual_hashes[file_path] != expected:
                raise ValueError(f"sha256 mismatch for {field} of {method!r}: {file_path}")
        hashes[method] = values
    missing = set(methods) - hashes.keys()
    if missing:
        raise ValueError(f"Missing hashes for methods: {', '.join(sorted(missing))}")
    if len({hashes[method] for method in methods}) != 1:
        raise ValueError("FASTQ/reference sha256 differ between methods")


def ratio(numerator, denominator):
    return f"{numerator / denominator:.6f}" if denominator else ""


def synthetic_tables(methods, truth):
    metrics = []
    depth_rows = []
    off_target_pairs = sum(off_target for _, off_target in truth.values())
    for method, calls in methods.items():
        depth = Counter(status if status != "classified" else RANKS[max(
            RANKS.index(rank) for rank in levels)]
            for levels, status in calls.values())
        for label in ("no_hit", "unclassified", *RANKS):
            depth_rows.append((method, label, depth[label], ratio(depth[label], len(calls))))
        for rank in EVALUATED_RANKS:
            tp = fp = fn = wrong_species = off_target_calls = unscored = 0
            for read_id, (expected, off_target) in truth.items():
                called = calls[read_id][0]
                has_call = rank in called
                if off_target:
                    if has_call:
                        fp += 1
                        off_target_calls += 1
                elif rank in expected:
                    if has_call and taxon_key(called, rank) == taxon_key(expected, rank):
                        tp += 1
                    else:
                        fn += 1
                        if has_call:
                            fp += 1
                            if rank == "species":
                                wrong_species += 1
                elif has_call:
                    if rank == "species":
                        fp += 1
                        wrong_species += 1
                    else:
                        unscored += 1
            metrics.append((method, rank, tp, fp, fn, ratio(tp, tp + fp),
                            ratio(tp, tp + fn), wrong_species, unscored,
                            off_target_calls, off_target_pairs,
                            ratio(off_target_calls * 1_000_000, off_target_pairs)))
    return metrics, depth_rows


def normalized(name):
    return " ".join(name.replace("_", " ").split()).casefold()


def mock_expected(path, run_accession):
    rows = [row for row in read_tsv(path, ("run_accession", "expected_species",
            "source_families", "coverage_status", "exact_species_reference_records",
            "genus_reference_records", "family_reference_records"))
            if row["run_accession"] == run_accession]
    if not rows:
        raise ValueError(f"No expected taxa for run {run_accession!r}")
    seen = set()
    expected = {rank: set() for rank in MOCK_RANKS}
    for row in rows:
        species = normalized(row["expected_species"])
        if not species or species in seen:
            raise ValueError(f"Empty or duplicate expected species for run {run_accession!r}")
        seen.add(species)
        expected["species"].add(species)
        expected["genus"].add(species.split()[0])
        expected["family"].update(normalized(name) for name in
                                  row["source_families"].split(",") if name.strip())
    return rows, expected


def mock_tables(methods, members, expected, min_fragments):
    metrics = []
    detection = []
    unexpected = []
    for method, calls in methods.items():
        support = {rank: Counter() for rank in MOCK_RANKS}
        labels = {rank: {} for rank in MOCK_RANKS}
        for levels, _ in calls.values():
            for rank in MOCK_RANKS:
                if rank not in levels:
                    continue
                name = levels[rank]
                key = normalized(name)
                support[rank][key] += 1
                labels[rank].setdefault(key, name)
        for rank in MOCK_RANKS:
            found = {key for key, count in support[rank].items()
                     if count >= min_fragments}
            metrics.append((method, rank, min_fragments, len(expected[rank]),
                            len(found & expected[rank]), len(found - expected[rank])))
            for key, count in sorted(support[rank].items()):
                if key not in expected[rank]:
                    unexpected.append((method, rank, labels[rank][key], count,
                                       "yes" if count >= min_fragments else "no"))
        for member in members:
            species = member["expected_species"]
            count = support["species"][normalized(species)]
            detection.append((method, member["run_accession"], species,
                              member["coverage_status"],
                              member["exact_species_reference_records"], count,
                              "yes" if count >= min_fragments else "no"))
    return metrics, detection, unexpected


def write_tsv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calls", action="append", required=True, metavar="METHOD=PATH")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--truth", type=Path)
    source.add_argument("--expected-taxa", type=Path)
    parser.add_argument("--run-accession")
    parser.add_argument("--hashes-manifest", type=Path)
    parser.add_argument("--min-fragments", type=int, default=1)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.min_fragments < 1:
        parser.error("--min-fragments must be positive")
    if args.expected_taxa and not args.run_accession:
        parser.error("--run-accession is required with --expected-taxa")
    methods = {}
    try:
        for item in args.calls:
            if "=" not in item:
                raise ValueError(f"Invalid --calls {item!r}; expected METHOD=PATH")
            name, path = item.split("=", 1)
            if not name or not path or name in methods:
                raise ValueError(f"Empty or duplicate method in --calls {item!r}")
            methods[name] = load_calls(Path(path))
        if args.hashes_manifest:
            check_hashes(args.hashes_manifest, methods)
        first_ids = set(next(iter(methods.values())))
        if not first_ids or any(set(calls) != first_ids for calls in methods.values()):
            raise ValueError("Methods have empty or differing ReadID sets")
        if args.truth:
            truth = load_truth(args.truth)
            if set(truth) != first_ids:
                raise ValueError("Truth and calls have differing ReadID sets")
            outputs = (
                ("rank_metrics.tsv", ("method", "rank", "TP", "FP", "FN", "precision",
                 "recall", "wrong_species_overcalls", "unscored_calls",
                 "off_target_false_calls", "off_target_pairs",
                 "false_calls_per_million_off_target_pairs")),
                ("resolution_depth.tsv", ("method", "depth", "pairs", "fraction")),
            )
            tables = synthetic_tables(methods, truth)
        else:
            members, expected = mock_expected(args.expected_taxa, args.run_accession)
            outputs = (
                ("mock_metrics.tsv", ("method", "rank", "min_fragments", "expected_taxa",
                 "detected_taxa", "unexpected_taxa")),
                ("expected_detection.tsv", ("method", "run_accession", "expected_species",
                 "coverage_status", "exact_species_reference_records", "fragments", "detected")),
                ("unexpected_calls.tsv", ("method", "rank", "taxon", "fragments",
                 "passes_threshold")),
            )
            tables = mock_tables(methods, members, expected, args.min_fragments)
        args.out.mkdir(parents=True, exist_ok=True)
        for (filename, fields), rows in zip(outputs, tables):
            write_tsv(args.out / filename, fields, rows)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
