#!/usr/bin/env python3
"""Include every classified pilot pair and sample no-hit pairs for BLASTn.

Each output file contains one ReadID per line for selecting both FASTQ mates.
Selection is independent of FASTQ and taxonomy row order.
"""

import argparse
import csv
import gzip
import hashlib
import json
import random
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fastq_ids(path):
    ids = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        while header := handle.readline():
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not header.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError(f"Malformed FASTQ record in {path}")
            if len(sequence.strip()) != len(quality.strip()):
                raise ValueError(f"Sequence/quality length mismatch in {path}")
            tokens = header[1:].split()
            if not tokens:
                raise ValueError(f"Empty read ID in {path}")
            read_id = tokens[0].removesuffix("/1")
            if not read_id:
                raise ValueError(f"Empty read ID in {path}")
            if read_id in ids:
                raise ValueError(f"Duplicate read-pair ID: {read_id}")
            ids.add(read_id)
    return ids


def read_taxonomy(path, read_ids):
    counts = {name: 0 for name in ("hit_above_genus", "hit_genus", "hit_species")}
    hit_ids = set()
    seen = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not {"ReadID", "Taxonomy"} <= set(header):
            raise ValueError("Taxonomy TSV must have ReadID and Taxonomy columns")
        for line, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed taxonomy TSV row at line {line}")
            read_id = row["ReadID"]
            if not read_id:
                raise ValueError(f"Empty ReadID at line {line}")
            if read_id in seen:
                raise ValueError(f"Duplicate taxonomy ReadID: {read_id}")
            if read_id not in read_ids:
                raise ValueError(f"Taxonomy ReadID not present in R1: {read_id}")
            seen.add(read_id)
            taxonomy = row["Taxonomy"].strip()
            if not taxonomy:
                continue
            depth = 0
            for level in taxonomy.split(";"):
                if not level or level.endswith("__"):
                    break
                depth += 1
            stratum = ("hit_above_genus" if depth < 6 else
                       "hit_genus" if depth == 6 else "hit_species")
            counts[stratum] += 1
            hit_ids.add(read_id)
    return hit_ids, counts


def write_ids(path, ids):
    ordered = sorted(ids)
    path.write_text("\n".join(ordered) + ("\n" if ordered else ""), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastq-r1", required=True, type=Path, help="Fastp-retained R1 FASTQ")
    parser.add_argument("--taxonomy", required=True, type=Path,
                        help="BarcodeMapper per-fragment ReadID/Taxonomy TSV")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--no-hit-per-seed", required=True, type=int)
    parser.add_argument("--seeds", nargs="+", type=int, required=True)
    args = parser.parse_args(argv)
    if args.no_hit_per_seed < 0:
        parser.error("--no-hit-per-seed must be nonnegative")
    if len(args.seeds) < 3 or len(set(args.seeds)) != len(args.seeds):
        parser.error("--seeds must contain at least three distinct integers")

    fastq = args.fastq_r1.expanduser().resolve()
    taxonomy = args.taxonomy.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    try:
        read_ids = read_fastq_ids(fastq)
        hit_ids, counts = read_taxonomy(taxonomy, read_ids)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    no_hit_ids = sorted(read_ids - hit_ids)
    counts["no_hit"] = len(no_hit_ids)
    if args.no_hit_per_seed > len(no_hit_ids):
        parser.error(f"There are only {len(no_hit_ids)} no-hit read pairs; "
                     f"requested {args.no_hit_per_seed}")
    if out_dir.exists():
        parser.error(f"Output directory already exists: {out_dir}")

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir()
    subsets = []
    union_ids = set()
    for seed in args.seeds:
        selected = hit_ids | set(random.Random(seed).sample(no_hit_ids, args.no_hit_per_seed))
        union_ids.update(selected)
        path = out_dir / f"subset.seed-{seed}.ids.txt"
        write_ids(path, selected)
        subsets.append({
            "seed": seed,
            "path": path.name,
            "read_pairs": len(selected),
            "no_hit_drawn": args.no_hit_per_seed,
            "sha256": sha256(path),
        })
    union_path = out_dir / "union.ids.txt"
    write_ids(union_path, union_ids)
    manifest = {
        "fastq_r1_path": str(fastq),
        "fastq_r1_sha256": sha256(fastq),
        "taxonomy_path": str(taxonomy),
        "taxonomy_sha256": sha256(taxonomy),
        "seeds": args.seeds,
        "no_hit_per_seed": args.no_hit_per_seed,
        "strata": [{"name": name, "available_pairs": counts[name]} for name in
                   ("hit_above_genus", "hit_genus", "hit_species", "no_hit")],
        "subsets": subsets,
        "union": {
            "path": union_path.name,
            "read_pairs": len(union_ids),
            "sha256": sha256(union_path),
        },
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(subsets)} subsets and union of {len(union_ids)} read pairs to {out_dir}")


if __name__ == "__main__":
    main()
