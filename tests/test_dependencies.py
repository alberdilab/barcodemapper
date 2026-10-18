"""Offline checks for the conda dependency checker."""

import importlib.util
import io
import json
import urllib.error
from pathlib import Path

import pytest


@pytest.fixture
def checker():
    script = Path(__file__).resolve().parents[1] / "scripts" / "check_dependencies.py"
    spec = importlib.util.spec_from_file_location("check_dependencies", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def env_file(checker, monkeypatch, tmp_path):
    source = ("name: barcodemapper\n"
              "channels:\n"
              "    - conda-forge\n"
              "    - bioconda\n"
              "dependencies:\n"
              "    - python=3.12\n"
              "    - fastp=0.23.1\n"
              "    - bowtie2=2.5.4\n")
    (tmp_path / "barcodemapper_env.yml").write_text(source, encoding="utf-8")
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    return source


def test_parses_environment_without_yaml(checker):
    text = (checker.ROOT / "barcodemapper_env.yml").read_text(encoding="utf-8")
    channels, pins = checker.parse_environment(text)
    assert channels == ["conda-forge", "bioconda"]
    assert pins["python"] == "3.12"
    assert pins["snakemake-executor-plugin-slurm"] == "1.1.0"
    assert pins["pyyaml"] == "6.0.2"
    assert len(pins) == 12


@pytest.mark.parametrize(("older", "newer"), [
    ("0.23.1", "0.24.0"),
    ("2.5.4", "2.10.0"),
    ("1.9", "1.85"),
])
def test_compares_dotted_versions_numerically(checker, older, newer):
    assert checker.compare_versions(newer, older) > 0
    assert checker.compare_versions(older, newer) < 0


def test_version_comparison_falls_back_to_strings(checker):
    assert checker.compare_versions("1.0b", "1.0a") > 0
    assert checker.compare_versions("1.0.0", "1.0") == 0


def test_fetch_latest_retries_with_timeout_and_user_agent(checker, monkeypatch):
    requests = []

    def open_url(request, timeout):
        requests.append((request, timeout))
        if len(requests) == 1:
            raise urllib.error.URLError("temporary failure")
        return io.BytesIO(b'{"latest_version": "0.24.0"}')

    monkeypatch.setattr(checker.urllib.request, "urlopen", open_url)
    assert checker.fetch_latest("bioconda", "fastp") == "0.24.0"
    assert len(requests) == 2
    assert requests[0][0].full_url == "https://api.anaconda.org/package/bioconda/fastp"
    assert requests[0][0].get_header("User-agent")
    assert requests[0][1] == 10


def test_exact_candidate_keeps_python_minor_and_pins_latest(checker, env_file, monkeypatch, tmp_path, capsys):
    releases = {
        ("conda-forge", "python"): "3.13.4",
        ("conda-forge", "fastp"): "0.24.0",
        ("bioconda", "fastp"): "0.23.2",
        ("conda-forge", "bowtie2"): "2.10.0",
        ("bioconda", "bowtie2"): "2.9.0",
    }
    monkeypatch.setattr(checker, "fetch_latest", lambda channel, name: releases.get((channel, name)))
    candidate = tmp_path / "candidate.yml"

    assert checker.main(["--write-candidate", str(candidate), "--exact"]) == 0
    output = candidate.read_text(encoding="utf-8")
    assert output == env_file.replace("name: barcodemapper\n", "name: barcodemapper-candidate\n") \
        .replace("fastp=0.23.1", "fastp=0.24.0").replace("bowtie2=2.5.4", "bowtie2=2.10.0")
    assert "newer available" in capsys.readouterr().out

    assert checker.main(["--write-candidate", str(candidate), "--exact", "--bump-python"]) == 0
    assert "python=3.13.4" in candidate.read_text(encoding="utf-8")


def test_default_candidate_leaves_versions_to_the_solver(checker, env_file, monkeypatch, tmp_path):
    releases = {("conda-forge", "python"): "3.13.4", ("conda-forge", "fastp"): "0.24.0", ("conda-forge", "bowtie2"): "2.5.4"}
    monkeypatch.setattr(checker, "fetch_latest", lambda channel, name: releases.get((channel, name)))
    candidate = tmp_path / "candidate.yml"
    assert checker.main(["--write-candidate", str(candidate)]) == 0
    dependencies = candidate.read_text(encoding="utf-8").split("dependencies:")[1]
    assert "- fastp\n" in dependencies and "fastp=" not in dependencies
    assert "python=3.12" in dependencies  # the Python minor stays pinned


def test_json_output_and_exit_codes(checker, env_file, monkeypatch, capsys):
    releases = {"python": "3.12", "fastp": "0.24.0", "bowtie2": "2.5.4"}
    monkeypatch.setattr(checker, "fetch_latest", lambda channel, name: releases[name] if channel == "conda-forge" else None)

    assert checker.main(["--json"]) == 0
    results = json.loads(capsys.readouterr().out)
    assert [row["package"] for row in results["outdated"]] == ["fastp"]
    assert [row["package"] for row in results["current"]] == ["python", "bowtie2"]
    assert results["failed"] == []
    assert checker.main(["--json", "--fail-if-outdated"]) == 10
    capsys.readouterr()

    def missing(channel, name):
        if name == "bowtie2":
            raise OSError("API unavailable")
        return releases[name] if channel == "conda-forge" else None

    monkeypatch.setattr(checker, "fetch_latest", missing)
    assert checker.main(["--json", "--fail-if-outdated"]) == 2
    results = json.loads(capsys.readouterr().out)
    assert results["failed"] == [{"package": "bowtie2", "pinned": "2.5.4",
                                  "latest": None, "status": "lookup failed"}]


def test_python_is_held_unless_bumped(checker, env_file, monkeypatch, capsys):
    releases = {"python": "3.13.4", "fastp": "0.23.1", "bowtie2": "2.5.4"}
    monkeypatch.setattr(checker, "fetch_latest", lambda channel, name: releases.get(name))
    assert checker.main(["--fail-if-outdated"]) == 0
    assert "held" in capsys.readouterr().out
    assert checker.main(["--fail-if-outdated", "--bump-python"]) == 10
