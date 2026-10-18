#!/usr/bin/env python3
"""Sample equal numbers of read pairs per truth-table stratum for BLASTn.

Example:
  make_blast_subset.py --truth truth.tsv --out-dir blast-subsets \
      --strata Marker Stratum --per-stratum 100 --seeds 1337 2027 4099

Each output file contains one ReadID per line for selecting both FASTQ mates.
The input must have exactly one row per read pair and enough rows in every
stratum. Selection is independent of input row order.
"""

import argparse
import csv
import hashlib
import json
import random
from collections import defaultdict
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_strata(path, id_column, columns):
    groups = defaultdict(list)
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)):
            raise ValueError("Truth TSV has duplicate column names")
        missing = [name for name in (id_column, *columns) if name not in header]
        if missing:
            raise ValueError(f"Truth TSV is missing columns: {', '.join(missing)}")
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed truth TSV row at line {line}")
            read_id = row[id_column]
            if not read_id:
                raise ValueError(f"Empty {id_column} at line {line}")
            if read_id in seen:
                raise ValueError(f"Duplicate read-pair ID: {read_id}")
            key = tuple(row[name] for name in columns)
            if any(not value for value in key):
                raise ValueError(f"Empty stratum value at line {line}")
            seen.add(read_id)
            groups[key].append(read_id)
    if not groups:
        raise ValueError("Truth TSV has no read pairs")
    return groups


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--truth", required=True, type=Path, help="One row per read pair, TSV")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--id-column", default="ReadID")
    parser.add_argument("--strata", nargs="+", required=True, metavar="COLUMN")
    parser.add_argument("--per-stratum", type=int, required=True)
    parser.add_argument("--seeds", nargs="+", type=int, required=True)
    args = parser.parse_args(argv)
    if args.per_stratum < 1:
        parser.error("--per-stratum must be positive")
    if len(args.seeds) < 3 or len(set(args.seeds)) != len(args.seeds):
        parser.error("--seeds must contain at least three distinct integers")
    if len(set(args.strata)) != len(args.strata):
        parser.error("--strata columns must be distinct")

    truth = args.truth.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    try:
        groups = read_strata(truth, args.id_column, args.strata)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    for key, ids in sorted(groups.items()):
        if len(ids) < args.per_stratum:
            parser.error(f"Stratum {key!r} has only {len(ids)} read pairs; "
                         f"requested {args.per_stratum}")
    if out_dir.exists():
        parser.error(f"Output directory already exists: {out_dir}")

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir()
    subsets = []
    for seed in args.seeds:
        rng = random.Random(seed)
        selected = []
        for key in sorted(groups):
            selected.extend(rng.sample(sorted(groups[key]), args.per_stratum))
        path = out_dir / f"subset.seed-{seed}.ids.txt"
        path.write_text("\n".join(sorted(selected)) + "\n", encoding="utf-8")
        subsets.append({
            "seed": seed,
            "path": path.name,
            "read_pairs": len(selected),
            "sha256": sha256(path),
        })
    manifest = {
        "input_path": str(truth),
        "input_sha256": sha256(truth),
        "id_column": args.id_column,
        "strata_columns": args.strata,
        "per_stratum": args.per_stratum,
        "strata": [{"values": list(key), "available_pairs": len(groups[key])}
                   for key in sorted(groups)],
        "subsets": subsets,
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(subsets)} subsets of {len(selected)} read pairs to {out_dir}")


if __name__ == "__main__":
    main()
