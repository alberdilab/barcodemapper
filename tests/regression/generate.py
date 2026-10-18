"""Build a small, byte-for-byte reproducible paired-read regression fixture."""

import argparse
import gzip
import random
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent / "data"
ADAPTER = "CTGTCTCTTATACACATCT"
ANIMAL = "k__Animalia;p__Arthropoda;c__Insecta;o__Lepidoptera"
FUNGUS = "k__Fungi;p__Ascomycota;c__Leotiomycetes;o__Helotiales"


def _sequence(seed: int, length: int = 220) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _mutate(sequence: str, positions: tuple[int, ...]) -> str:
    bases = list(sequence)
    for position in positions:
        bases[position] = "C" if bases[position] == "A" else "A"
    return "".join(bases)


def _reads(sequence: str) -> tuple[str, str]:
    return sequence[:125], _reverse_complement(sequence[-125:])


def _write_sample(directory: Path, sample: str, pairs: list[tuple[str, str]]) -> None:
    for mate in (1, 2):
        with (directory / f"{sample}_{mate}.fq.gz").open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=9, mtime=0, fileobj=raw) as handle:
                for index, pair in enumerate(pairs, 1):
                    read = pair[mate - 1]
                    handle.write(
                        f"@{sample}_{index:03d}\n{read}\n+\n{'I' * len(read)}\n".encode("ascii")
                    )


def generate_fixture(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    reads_dir = output / "reads"
    reads_dir.mkdir(exist_ok=True)

    references = [
        ("COI01", "COI-5P", f"{ANIMAL};f__Hesperiidae;g__Solitaria;s__alpha", 101),
        ("COI02", "COI-5P", f"{ANIMAL};f__Nymphalidae;g__Ambigua;s__one", 102),
        ("COI03", "COI-5P", f"{ANIMAL};f__Nymphalidae;g__Ambigua;s__two", 102),
        ("COI04", "COI-5P", f"{ANIMAL};f__Lycaenidae;g__Gemina;s__one", 103),
        ("COI05", "COI-5P", f"{ANIMAL};f__Papilionidae;g__Gemina;s__two", 103),
        ("COI06", "COI-5P", f"{ANIMAL};f__Pieridae;g__Mutata;s__alpha", 104),
        ("ITS01", "ITS", f"{FUNGUS};f__Sclerotiniaceae;g__Solifungus;s__alpha", 201),
        ("ITS02", "ITS", f"{FUNGUS};f__Helotiaceae;g__Ambifungus;s__one", 202),
        ("ITS03", "ITS", f"{FUNGUS};f__Helotiaceae;g__Ambifungus;s__two", 202),
        ("ITS04", "ITS", f"{FUNGUS};f__Hyaloscyphaceae;g__Gemifungus;s__one", 203),
        ("ITS05", "ITS", f"{FUNGUS};f__Dermateaceae;g__Gemifungus;s__two", 203),
        ("ITS06", "ITS", f"{FUNGUS};f__Rutstroemiaceae;g__Mutafungus;s__alpha", 204),
    ]
    sequences = {accession: _sequence(seed) for accession, _, _, seed in references}
    (output / "references.fa").write_text(
        "".join(
            f">{accession}|{marker}|{taxonomy}\n{sequences[accession]}\n"
            for accession, marker, taxonomy, _ in references
        ),
        encoding="ascii",
    )

    # Each mate has three edits (within 3% of 125 bp) or four (just beyond it).
    within = (20, 50, 80, 140, 170, 200)
    beyond = (20, 40, 60, 80, 140, 160, 180, 200)

    def mapped(accession: str, count: int) -> list[tuple[str, str]]:
        return [_reads(sequences[accession])] * count

    def edited(accession: str, positions: tuple[int, ...], count: int) -> list[tuple[str, str]]:
        return [_reads(_mutate(sequences[accession], positions))] * count

    def noise(count: int, seed: int) -> list[tuple[str, str]]:
        return [_reads(_sequence(seed + index)) for index in range(count)]

    def filtered(count: int, seed: int) -> list[tuple[str, str]]:
        low_complexity = [("A" * 125, "T" * 125)] * count
        adaptered = [
            (_sequence(seed + index, 90) + ADAPTER, _sequence(seed + 100 + index, 90) + ADAPTER)
            for index in range(count)
        ]
        return low_complexity + adaptered

    samples = {
        "animal": (mapped("COI01", 8) + mapped("COI02", 6) + mapped("COI04", 4)
                   + edited("COI06", within, 3) + edited("COI06", beyond, 3)
                   + noise(3, 300) + filtered(3, 400)),
        "fungal": (mapped("ITS01", 8) + mapped("ITS02", 6) + mapped("ITS04", 4)
                   + edited("ITS06", within, 3) + edited("ITS06", beyond, 3)
                   + noise(3, 500) + filtered(3, 600)),
        "mixed": (mapped("COI01", 3) + mapped("ITS01", 3) + mapped("COI02", 2)
                  + mapped("ITS02", 2) + mapped("COI04", 2) + mapped("ITS04", 2)
                  + edited("COI06", within, 2) + edited("ITS06", within, 2)
                  + noise(4, 700) + filtered(3, 800)),
        "zero": noise(12, 900) + filtered(6, 1000),
    }
    for sample, pairs in samples.items():
        _write_sample(reads_dir, sample, pairs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DATA_DIR, help="Fixture data directory")
    args = parser.parse_args()
    generate_fixture(args.output)


if __name__ == "__main__":
    main()
