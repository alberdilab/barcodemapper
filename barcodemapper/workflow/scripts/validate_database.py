"""Validation of the reference FASTA consumed by Bowtie2 and taxonomy assignment."""

import argparse
import re
from pathlib import Path


RANKS = {rank: index for index, rank in enumerate("kpcofgs")}
RANK_PATTERN = re.compile(r"^([kpcofgs])__(.*)$")


def validate_database(path: str | Path) -> dict[str, int]:
    """Check reference identifiers and return the number of records per marker.

    Bowtie2 uses only the first whitespace-delimited token of a FASTA header as
    its reference name, so headers with whitespace cannot safely carry taxonomy.
    """
    counts = {"COI-5P": 0, "ITS": 0}
    records = 0
    sequence_bases = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if records and not sequence_bases:
                    raise ValueError(f"Reference record before line {line_number} has no sequence")
                header = line[1:]
                if any(char.isspace() for char in header):
                    raise ValueError(f"Whitespace in reference header at line {line_number}")
                parts = header.split("|")
                if len(parts) != 3 or not parts[0]:
                    raise ValueError(f"Expected accession|marker|taxonomy at line {line_number}")
                marker = parts[1]
                if marker not in counts:
                    raise ValueError(f"Unsupported marker {marker!r} at line {line_number}")
                previous_rank = -1
                for index, level in enumerate(parts[2].split(";")):
                    match = RANK_PATTERN.fullmatch(level)
                    if match is None:
                        raise ValueError(f"Invalid taxonomy rank at line {line_number}: {level!r}")
                    rank, name = match.groups()
                    rank_index = RANKS[rank]
                    if (index == 0 and (rank != "k" or not name)) or rank_index <= previous_rank:
                        raise ValueError(f"Taxonomy ranks are missing a kingdom or out of order at line {line_number}")
                    previous_rank = rank_index
                counts[marker] += 1
                records += 1
                sequence_bases = 0
            elif not records:
                raise ValueError(f"Sequence before the first reference header at line {line_number}")
            else:
                sequence_bases += len(line)
    if not records:
        raise ValueError("Reference FASTA contains no records")
    if not sequence_bases:
        raise ValueError("Last reference record has no sequence")
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a BarcodeMapper reference FASTA")
    parser.add_argument("fasta")
    args = parser.parse_args()
    validate_database(args.fasta)
