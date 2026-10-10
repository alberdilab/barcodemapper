"""Generate a tiny deterministic paired-read dataset for workflow checks."""

import gzip
import random
from pathlib import Path


TAXON_PREFIX = "k__Animalia;p__Chordata;c__Mammalia;o__Primates;f__Hominidae;g__Synthus"


def _sequence(seed: int, length: int = 220) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_pair(directory: Path, sample: str, sequence: str) -> None:
    read1 = sequence[:125]
    read2 = _reverse_complement(sequence[-125:])
    for mate, read in ((1, read1), (2, read2)):
        with gzip.open(directory / f"{sample}_{mate}.fq.gz", "wt") as handle:
            handle.write(f"@{sample}\n{read}\n+\n{'I' * len(read)}\n")


def make_fixture(root: Path) -> tuple[Path, Path]:
    reads = root / "reads"
    reads.mkdir(parents=True)
    target = _sequence(712)
    database = root / "synthetic.fa"
    database.write_text(
        f">ref_A|COI-5P|{TAXON_PREFIX};s__Synthus_one\n{target}\n"
        f">ref_B|COI-5P|{TAXON_PREFIX};s__Synthus_two\n{target}\n",
        encoding="utf-8",
    )
    _write_pair(reads, "positive", target)
    _write_pair(reads, "negative", _sequence(984))
    return database, reads
