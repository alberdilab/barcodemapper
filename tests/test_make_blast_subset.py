"""Offline checks for reproducible BLAST read-pair subsets."""

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/mjolnir/make_blast_subset.py"


def write_truth(path):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("ReadID", "Marker", "Stratum", "Taxonomy"))
        for marker in ("COI", "ITS"):
            for stratum in ("present", "absent", "ambiguous", "off_target"):
                for number in range(4):
                    writer.writerow((f"{marker}-{stratum}-{number}", marker, stratum, ""))


def run_sampler(truth, output, per_stratum=2, seeds=(11, 23, 37)):
    return subprocess.run([
        sys.executable, str(SCRIPT), "--truth", str(truth), "--out-dir", str(output),
        "--strata", "Marker", "Stratum", "--per-stratum", str(per_stratum),
        "--seeds", *(str(seed) for seed in seeds),
    ], capture_output=True, text=True)


def test_same_input_and_seeds_produce_identical_files(tmp_path):
    truth = tmp_path / "truth.tsv"
    write_truth(truth)
    first = tmp_path / "first"
    second = tmp_path / "second"
    assert run_sampler(truth, first).returncode == 0
    assert run_sampler(truth, second).returncode == 0

    names = sorted(path.name for path in first.iterdir())
    assert names == ["manifest.json", "subset.seed-11.ids.txt", "subset.seed-23.ids.txt",
                     "subset.seed-37.ids.txt"]
    for name in names:
        assert (first / name).read_bytes() == (second / name).read_bytes()

    manifest = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["input_sha256"] == hashlib.sha256(truth.read_bytes()).hexdigest()
    for subset in manifest["subsets"]:
        data = (first / subset["path"]).read_bytes()
        assert subset["sha256"] == hashlib.sha256(data).hexdigest()


def test_each_seed_covers_every_marker_and_challenge_stratum(tmp_path):
    truth = tmp_path / "truth.tsv"
    write_truth(truth)
    output = tmp_path / "subsets"
    assert run_sampler(truth, output).returncode == 0
    for path in output.glob("*.ids.txt"):
        ids = path.read_text(encoding="utf-8").splitlines()
        assert len(ids) == len(set(ids)) == 16
        assert ids == sorted(ids)
        for marker in ("COI", "ITS"):
            for stratum in ("present", "absent", "ambiguous", "off_target"):
                assert sum(read_id.startswith(f"{marker}-{stratum}-") for read_id in ids) == 2


def test_rejects_underfilled_strata_without_writing_subsets(tmp_path):
    truth = tmp_path / "truth.tsv"
    write_truth(truth)
    output = tmp_path / "subsets"
    result = run_sampler(truth, output, per_stratum=5)
    assert result.returncode != 0
    assert "only 4" in result.stderr
    assert not output.exists()


def test_requires_three_distinct_seeds(tmp_path):
    truth = tmp_path / "truth.tsv"
    write_truth(truth)
    output = tmp_path / "subsets"
    result = run_sampler(truth, output, seeds=(11, 23))
    assert result.returncode != 0
    assert "at least three distinct" in result.stderr
    assert not output.exists()
