"""Regression fixture and real-workflow golden checks."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "regression" / "data"
GENERATOR = ROOT / "tests" / "regression" / "generate.py"
CHECK = ROOT / "scripts" / "regression_check.py"


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
def test_regression_taxonomy_matches_golden():
    tools = ("snakemake", "fastp", "bowtie2", "bowtie2-build", "samtools")
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        pytest.skip("workflow executables unavailable: " + ", ".join(missing))
    result = subprocess.run(
        [sys.executable, str(CHECK)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
