"""Regression fixture and real-workflow golden checks."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from barcodemapper import cli


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "regression" / "data"
GENERATOR = ROOT / "tests" / "regression" / "generate.py"
GOLDEN = ROOT / "tests" / "regression" / "expected_taxonomy.tsv"


def _fixture_bytes(directory: Path) -> dict[str, bytes]:
    return {
        str(path.relative_to(directory)): path.read_bytes()
        for path in directory.rglob("*") if path.is_file()
    }


def test_regression_fixture_is_deterministic(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    for directory in (first, second):
        subprocess.run([sys.executable, str(GENERATOR), "--output", str(directory)], check=True)
    assert _fixture_bytes(first) == _fixture_bytes(second) == _fixture_bytes(DATA)


@pytest.mark.integration
def test_regression_taxonomy_matches_golden(tmp_path):
    tools = ("snakemake", "fastp", "bowtie2", "bowtie2-build", "samtools")
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        pytest.skip("workflow executables unavailable: " + ", ".join(missing))
    output = tmp_path / "taxonomy.tsv"
    cli.main([
        "-i", str(DATA / "reads"), "-d", str(DATA / "references.fa"),
        "-o", str(output), "-t", str(tmp_path / "work"), "--standard",
    ])
    actual = output.read_text().splitlines()
    expected = GOLDEN.read_text().splitlines()
    assert actual[0] == expected[0]
    assert sorted(actual[1:]) == sorted(expected[1:])
