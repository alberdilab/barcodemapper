#!/usr/bin/env python3
"""Compare per-read taxonomy assignments of capped (-k N) runs with an all-hits (-a) run.

Assignments are reduced to their resolved prefix (highest contiguous named ranks) first,
so they match what the pipeline reports in its final table. Reads only local files. Per sample, pass the -a baseline taxonomy table and one
table per cap; every table is the pipeline's `taxonomy/<sample>.tsv` (ReadID, Taxonomy).
Optional --resources adds wall time and memory per run from a TSV with columns
sample, label, elapsed_seconds, max_rss_mb.

Example:
  compare_max_hits.py --sample animal --baseline a=ERR6001161.taxonomy.tsv \
      --run k5=k5.taxonomy.tsv --run k50=k50.taxonomy.tsv --out comparison.tsv
"""
import argparse
import csv
import sys

RANKS = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order",
         "f": "family", "g": "genus", "s": "species"}


def resolved_prefix(taxonomy):
    """Highest contiguous named ranks, as the pipeline's aggregation step keeps them."""
    levels = []
    for level in (taxonomy or "").split(";"):
        if not level or level.endswith("__"):
            break
        levels.append(level)
    return ";".join(levels)


def read_assignments(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        return {row[0]: resolved_prefix(row[1]) for row in reader if len(row) >= 2}


def terminal_rank(taxonomy):
    return RANKS.get(taxonomy.rsplit(";", 1)[-1][:1], "unknown") if taxonomy else "none"


def relation(baseline, test):
    """How a capped assignment relates to the all-hits one for the same read."""
    if test == baseline:
        return "identical"
    if not test:
        return "less_specific"
    if baseline.startswith(test + ";"):
        return "less_specific"
    if test.startswith(baseline + ";") or not baseline:
        return "more_specific"
    return "different_branch"


def compare(baseline, test):
    counts = {"reads_baseline": len(baseline), "reads_run": len(test),
              "lost": 0, "gained": 0, "identical": 0, "less_specific": 0,
              "more_specific": 0, "different_branch": 0}
    for read, taxonomy in baseline.items():
        if read not in test:
            counts["lost"] += 1
        else:
            counts[relation(taxonomy, test[read])] += 1
    counts["gained"] = sum(read not in baseline for read in test)
    for label, table in (("base", baseline), ("run", test)):
        for rank in ("species", "genus", "family"):
            counts[f"{label}_resolved_{rank}"] = sum(terminal_rank(t) == rank for t in table.values())
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--baseline", required=True, help="LABEL=PATH of the -a taxonomy table")
    parser.add_argument("--run", action="append", default=[], help="LABEL=PATH of a capped run (repeatable)")
    parser.add_argument("--resources", help="Optional TSV: sample, label, elapsed_seconds, max_rss_mb")
    parser.add_argument("--out", help="Write the table here as TSV (default: stdout)")
    args = parser.parse_args(argv)

    def split(spec):
        label, _, path = spec.partition("=")
        if not path:
            parser.error(f"expected LABEL=PATH, got {spec!r}")
        return label, path

    resources = {}
    if args.resources:
        with open(args.resources, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                resources[(row["sample"], row["label"])] = row
    base_label, base_path = split(args.baseline)
    baseline = read_assignments(base_path)
    rows = []
    for spec in [args.baseline] + args.run:
        label, path = split(spec)
        row = {"sample": args.sample, "run": label, **compare(baseline, read_assignments(path))}
        extra = resources.get((args.sample, label), {})
        row["elapsed_seconds"] = extra.get("elapsed_seconds", "")
        row["max_rss_mb"] = extra.get("max_rss_mb", "")
        rows.append(row)
    handle = open(args.out, "w", newline="", encoding="utf-8") if args.out else sys.stdout
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    if args.out:
        handle.close()


if __name__ == "__main__":
    main()
