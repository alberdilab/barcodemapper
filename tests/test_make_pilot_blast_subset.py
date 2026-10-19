"""Offline checks for pilot BLAST read-pair subsets."""

import gzip
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/mjolnir/make_pilot_blast_subset.py"
SEEDS = (11, 23, 37)


def write_fastq(path, headers):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8") as handle:
        for header in headers:
            handle.write(f"@{header}\nACGT\n+\n!!!!\n")


def run_sampler(fastq, taxonomy, output, no_hit_per_seed=3, seeds=SEEDS):
    return subprocess.run([
        sys.executable, str(SCRIPT), "--fastq-r1", str(fastq),
        "--taxonomy", str(taxonomy), "--out-dir", str(output),
        "--no-hit-per-seed", str(no_hit_per_seed),
        "--seeds", *(str(seed) for seed in seeds),
    ], capture_output=True, text=True)


@pytest.mark.parametrize("suffix", (".fq", ".fq.gz"))
def test_subsets_include_all_hits_and_sample_no_hits(tmp_path, suffix):
    fastq = tmp_path / f"reads{suffix}"
    headers = ["hit0", "hit5/1 lane", "hit6 lane:16752/1", "hit7"]
    headers += [f"miss{number}" for number in range(12)]
    write_fastq(fastq, headers)
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text(
        "ReadID\tTaxonomy\n"
        "hit0\tk__\n"
        "hit5\tk__A;p__B;c__C;o__D;f__E;g__;s__F\n"
        "hit6\tk__A;p__B;c__C;o__D;f__E;g__F;s__\n"
        "hit7\tk__A;p__B;c__C;o__D;f__E;g__F;s__G\n"
        "miss0\t\n", encoding="utf-8")
    first = tmp_path / "first"
    second = tmp_path / "second"
    assert run_sampler(fastq, taxonomy, first).returncode == 0
    assert run_sampler(fastq, taxonomy, second).returncode == 0

    names = sorted(path.name for path in first.iterdir())
    assert names == ["manifest.json", "subset.seed-11.ids.txt", "subset.seed-23.ids.txt",
                     "subset.seed-37.ids.txt", "union.ids.txt"]
    for name in names:
        assert (first / name).read_bytes() == (second / name).read_bytes()

    manifest = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["fastq_r1_path"] == str(fastq.resolve())
    assert manifest["fastq_r1_sha256"] == hashlib.sha256(fastq.read_bytes()).hexdigest()
    assert manifest["taxonomy_path"] == str(taxonomy.resolve())
    assert manifest["taxonomy_sha256"] == hashlib.sha256(taxonomy.read_bytes()).hexdigest()
    assert manifest["seeds"] == list(SEEDS)
    assert manifest["no_hit_per_seed"] == 3
    assert manifest["strata"] == [
        {"name": "hit_above_genus", "available_pairs": 2},
        {"name": "hit_genus", "available_pairs": 1},
        {"name": "hit_species", "available_pairs": 1},
        {"name": "no_hit", "available_pairs": 12},
    ]

    hit_ids = {"hit0", "hit5", "hit6", "hit7"}
    no_hit_ids = sorted(f"miss{number}" for number in range(12))
    selected_by_seed = []
    for subset in manifest["subsets"]:
        seed = subset["seed"]
        path = first / subset["path"]
        selected = path.read_text(encoding="utf-8").splitlines()
        assert selected == sorted(selected)
        assert set(selected) == hit_ids | set(random.Random(seed).sample(no_hit_ids, 3))
        assert subset["read_pairs"] == 7
        assert subset["no_hit_drawn"] == 3
        assert subset["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        selected_by_seed.append(set(selected))
    assert len({tuple(sorted(ids)) for ids in selected_by_seed}) > 1
    union = first / "union.ids.txt"
    assert union.read_text(encoding="utf-8").splitlines() == sorted(set.union(*selected_by_seed))
    assert manifest["union"] == {
        "path": union.name,
        "read_pairs": len(set.union(*selected_by_seed)),
        "sha256": hashlib.sha256(union.read_bytes()).hexdigest(),
    }


def test_zero_hits_are_allowed(tmp_path):
    fastq = tmp_path / "reads.fq"
    write_fastq(fastq, ["r1", "r2", "r3", "r4"])
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text("ReadID\tTaxonomy\nr1\t\n", encoding="utf-8")
    output = tmp_path / "subsets"
    assert run_sampler(fastq, taxonomy, output, no_hit_per_seed=2).returncode == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert [group["available_pairs"] for group in manifest["strata"]] == [0, 0, 0, 4]
    for subset in manifest["subsets"]:
        assert subset["read_pairs"] == 2


@pytest.mark.parametrize("case,expected", [
    ("too_many", "only 2 no-hit"),
    ("duplicate", "Duplicate read-pair ID"),
    ("absent", "not present in R1"),
    ("few_seeds", "at least three distinct"),
    ("repeated_seeds", "at least three distinct"),
    ("existing", "already exists"),
])
def test_rejects_invalid_inputs_without_writing_subsets(tmp_path, case, expected):
    fastq = tmp_path / "reads.fq"
    write_fastq(fastq, ["r1/1", "r2", "r3"] if case != "duplicate" else ["r1", "r1/1", "r3"])
    taxonomy = tmp_path / "taxonomy.tsv"
    taxonomy.write_text("ReadID\tTaxonomy\n" +
                        ("missing\tk__A\n" if case == "absent" else "r1\tk__A\n"),
                        encoding="utf-8")
    output = tmp_path / "subsets"
    if case == "existing":
        output.mkdir()
    n = 3 if case == "too_many" else 1
    seeds = (11, 23) if case == "few_seeds" else (11, 23, 11) if case == "repeated_seeds" else SEEDS
    result = run_sampler(fastq, taxonomy, output, no_hit_per_seed=n, seeds=seeds)
    assert result.returncode != 0
    assert expected in result.stderr
    assert not output.exists() or list(output.iterdir()) == []
