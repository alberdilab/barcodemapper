"""Offline checks for portable database metadata and attribution."""

from copy import deepcopy
from getpass import getuser
from pathlib import Path

import pytest

from barcodemapper.publish_meta import attribution_text, portable_manifest, zenodo_metadata


@pytest.fixture
def automatic_manifest():
    return {
        "created_utc": "2026-10-21T12:00:00+00:00",
        "barcodemapper_version": "1.2.4",
        "database": {"path": "/work/bmdb_20261021/bmdb_20261021.fa.gz", "sha256": "database-hash"},
        "parameters": {
            "bold_retain": "k__Animalia",
            "unite_retain": "k__Fungi",
            "tmp_dir": "/work/tmp",
            "barcodemapper_db": "/work/bmdb_20261021/bmdb_20261021.fa.gz",
            "output_files": {"species": "/work/results/species.tsv"},
        },
        "samples": {"sample": {"read1": ["/work/reads/sample_1.fastq.gz"]}},
        "levels": [{"name": "species", "output_file": "/work/results/species.tsv"}],
        "tools": {"bowtie2": "2.5.4"},
        "automatic_sources": True,
        "sources": {
            "bold": {"release": "BOLD_Public.09-Oct-2026", "doi": "https://doi.org/10.5883/DP-Latest",
                     "url": "https://bench.boldsystems.org/fasta", "sha256": "bold-hash"},
            "unite": {"release": "2025-02-19", "doi": "https://doi.org/10.15156/BIO/3301232",
                      "url": "https://s3.hpc.ut.ee/unite.tgz", "sha256": "unite-hash"},
        },
    }


@pytest.fixture
def manual_manifest(automatic_manifest):
    manifest = deepcopy(automatic_manifest)
    manifest.pop("automatic_sources")
    manifest["sources"] = {
        "bold": {"path": "/home/user/bold.fa", "sha256": "bold-hash"},
        "unite": {"path": "/home/user/unite.fa", "sha256": "unite-hash"},
    }
    return manifest


def _strings(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from _strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _strings(item)
    elif isinstance(value, str):
        yield value


@pytest.mark.parametrize("fixture_name", ["automatic_manifest", "manual_manifest"])
def test_portable_manifest_scrubs_paths_without_changing_input(request, fixture_name):
    manifest = request.getfixturevalue(fixture_name)
    original = deepcopy(manifest)

    portable = portable_manifest(manifest, record_doi="10.5281/zenodo.123", record_url="https://zenodo.org/records/123")

    assert manifest == original
    assert not any(value.startswith("/") for value in _strings(portable))
    assert portable["database"] == {"path": "bmdb_20261021.fa.gz", "sha256": "database-hash"}
    assert portable["sources"]["bold"]["sha256"] == "bold-hash"
    assert portable["sources"]["unite"]["sha256"] == "unite-hash"
    assert portable["parameters"]["bold_retain"] == "k__Animalia"
    assert portable["parameters"]["unite_retain"] == "k__Fungi"
    assert portable["origin"] == "zenodo"
    assert portable["zenodo"] == {"doi": "10.5281/zenodo.123", "url": "https://zenodo.org/records/123"}
    if fixture_name == "automatic_manifest":
        assert portable["sources"] == original["sources"]
        assert portable["automatic_sources"] is True
    else:
        assert not portable["sources"]["bold"].get("path")
        assert not portable["sources"]["unite"].get("path")


def test_zenodo_metadata_and_attribution_include_source_credit(automatic_manifest):
    creators = [{"name": "Alberdi, Antton"}]
    metadata = zenodo_metadata(automatic_manifest, version="2026-10-21", creators=creators)

    assert metadata["upload_type"] == "dataset"
    assert metadata["license"] == "cc-by-sa-4.0"
    assert metadata["version"] == "2026-10-21"
    assert metadata["creators"] == creators
    assert "modified" in metadata["description"]
    assert "MIT" in metadata["description"]
    assert "k__Animalia" in metadata["description"]
    assert "k__Fungi" in metadata["description"]
    for source in automatic_manifest["sources"].values():
        assert source["release"] in metadata["description"]
        assert source["doi"] in metadata["description"]
        assert any(item["relation"] == "isDerivedFrom" and
                   item["identifier"] == source["doi"].removeprefix("https://doi.org/") and
                   item["scheme"] == "doi" for item in metadata["related_identifiers"])
    assert any(item["relation"] == "isSupplementTo" and item["scheme"] == "url" for item in metadata["related_identifiers"])

    attribution = attribution_text(automatic_manifest)
    assert "https://creativecommons.org/licenses/by-sa/4.0/" in attribution
    assert "changes made" in attribution.lower()
    assert "Ratnasingham & Hebert 2007" in attribution
    assert "Ratnasingham & Hebert 2013" in attribution
    for source in automatic_manifest["sources"].values():
        assert source["doi"] in attribution
        assert source["sha256"] in attribution


def test_manual_manifest_uses_hashes_without_source_dois(manual_manifest):
    metadata = zenodo_metadata(manual_manifest, version="manual", creators=[])
    attribution = attribution_text(manual_manifest)

    assert "bold-hash" in metadata["description"]
    assert "unite-hash" in metadata["description"]
    assert not any(item["relation"] == "isDerivedFrom" for item in metadata["related_identifiers"])
    assert "bold-hash" in attribution
    assert "unite-hash" in attribution


def test_description_escapes_manifest_text(automatic_manifest):
    automatic_manifest["sources"]["bold"]["release"] = "release <script>"
    automatic_manifest["parameters"]["bold_retain"] = "A & B"

    description = zenodo_metadata(automatic_manifest, version="1", creators=[])["description"]

    assert "release &lt;script&gt;" in description
    assert "A &amp; B" in description
    assert "<script>" not in description


def test_citations_dates_licences_and_alias(automatic_manifest):
    description = zenodo_metadata(automatic_manifest, version="20261021", creators=[])["description"]
    attribution = attribution_text(automatic_manifest)
    for text in (description, attribution):
        assert "BOLD_Public.09-Oct-2026" in text
        assert "Accessed/built: 2026-10-21" in text
        assert "10.5883/DP-Latest is a weekly-updating alias" in text
        assert "BOLD v4: A Centralized Bioinformatics Platform" in text
        assert "Ratnasingham" in text and "2007" in text and "2013" in text
        assert "Abarenkov K, Zirk A, Piirmann T, Pöhönen R, Ivanov F, Nilsson RH, Kõljalg U (2025)" in text
        assert "UNITE general FASTA release for eukaryotes 2" in text
        assert "https://creativecommons.org/licenses/by-sa/4.0/" in text
        assert "https://creativecommons.org/licenses/by/4.0/" in text
        assert "UNITE-derived portion is also available from UNITE under CC BY 4.0" in text


def test_portable_manifest_scrubs_embedded_paths_and_user(automatic_manifest):
    user = getuser()
    automatic_manifest["tools"] = {"bowtie2": f"2.5.4 built by {user} at /projects/alberdilab/scratch/jpl786/build",
                                   "/Users/private/tool": "secret"}
    automatic_manifest["parameters"].update({"db_dir": "/projects/alberdilab/scratch/jpl786/db",
                                               "threads": 8, "comment": f"built in ~/work by {user}"})
    automatic_manifest["sources"]["bold"]["note"] = "from C:\\private\\database and /home/user/raw"
    automatic_manifest["home"] = str(Path.home())
    serialized = __import__("json").dumps(portable_manifest(automatic_manifest))
    for private in ("/projects/", "/home/", "/Users/", "C:\\private", "~/", user, str(Path.home())):
        assert private not in serialized
    portable = portable_manifest(automatic_manifest)
    assert portable["parameters"]["threads"] == 8
    assert "db_dir" not in portable["parameters"]
    assert "/Users/private/tool" not in portable["tools"]
