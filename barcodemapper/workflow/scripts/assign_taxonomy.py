"""Assign each read ID to the taxonomy shared by its accepted alignments."""

import csv
from itertools import groupby


def common_taxonomy(first: str, second: str) -> str:
    shared = []
    for left, right in zip(first.split(";"), second.split(";")):
        if left != right:
            break
        shared.append(left)
    return ";".join(shared)


def assign_file(input_file: str, output_file: str) -> None:
    # The input BAM is query-name collated, so every alignment for a fragment is
    # adjacent in the extracted TSV. Process one fragment at a time.
    with open(input_file, newline="", encoding="utf-8") as source, open(
        output_file, "w", newline="", encoding="utf-8"
    ) as destination:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(("ReadID", "Taxonomy"))
        for read_id, rows in groupby(reader, key=lambda row: row[0]):
            taxonomy = None
            for row in rows:
                if len(row) != 6:
                    raise ValueError(f"Invalid extracted alignment for read {read_id!r}")
                taxonomy = row[3] if taxonomy is None else common_taxonomy(taxonomy, row[3])
            writer.writerow((read_id, taxonomy))


assign_file(snakemake.input[0], snakemake.output[0])
