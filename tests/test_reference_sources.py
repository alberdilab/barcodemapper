"""Offline checks for release discovery and compressed source handling."""

import gzip
import io
import json
import tarfile
from email.message import Message

import pytest

from barcodemapper import reference_sources as sources


class Response(io.BytesIO):
    def __init__(self, data: bytes, content_type: str = "text/html"):
        super().__init__(data)
        self.headers = Message()
        self.headers["Content-Type"] = content_type


def test_discovers_current_all_eukaryotes_with_singletons(monkeypatch):
    bold = b'<button name="download-fasta-unauthenticated" data-package-id="BOLD_Public.09-Oct-2026"></button>'
    unite = b"""<b>List of all general FASTA releases</b><table>
      <tr><td>10.0</td><td>2025-02-19</td><td>All eukaryotes</td><td>20 802</td><td>147 735</td><td>Current</td><td><a href="https://doi.org/10.15156/BIO/3301231">DOI</a></td><td>Without singletons</td></tr>
      <tr><td>10.0</td><td>2025-02-19</td><td>All eukaryotes</td><td>20 802</td><td>245 787</td><td>Current</td><td><a href="https://doi.org/10.15156/BIO/3301232">DOI</a></td><td>With singletons</td></tr>
      </table>"""
    doi = json.dumps({"data": [{"attributes": {"media": [{"name": "sh_general_release_s_all_19.02.2025.tgz", "url": "https://s3.hpc.ut.ee/unite.tgz"}]}}]}).encode()

    def open_source(url, **kwargs):
        if url == sources.BOLD_LATEST:
            return Response(bold)
        if url == sources.UNITE_REPOSITORY:
            return Response(unite)
        assert "identifier=10.15156%2FBIO%2F3301232" in url
        return Response(doi, "application/json")

    monkeypatch.setattr(sources, "_open", open_source)
    releases = sources.discover_latest()
    assert releases["bold"].identifier == "BOLD_Public.09-Oct-2026"
    assert releases["unite"].doi == "https://doi.org/10.15156/BIO/3301232"
    assert releases["unite"].url == "https://s3.hpc.ut.ee/unite.tgz"


def test_downloads_and_extracts_only_standard_unite_fasta(monkeypatch, tmp_path):
    bold_fasta = b">A|COI-5P|x|Animalia,Arthropoda\nACGT\n"
    unite_fasta = b">record|U1|x|x|k__Fungi\nTGCA\n"
    archive_bytes = io.BytesIO()
    with tarfile.open(fileobj=archive_bytes, mode="w:gz") as archive:
        for name, contents in (
            ("sh_general_release_dynamic_s_all_19.02.2025.fasta", unite_fasta),
            ("sh_general_release_dynamic_s_all_19.02.2025_dev.fasta", b"wrong"),
        ):
            member = tarfile.TarInfo(name)
            member.size = len(contents)
            archive.addfile(member, io.BytesIO(contents))

    def open_source(url, **kwargs):
        if "uid=" not in url and url.startswith(sources.BOLD_FASTA):
            return Response(b'"generateduid"', "application/json")
        if "uid=" in url:
            return Response(gzip.compress(bold_fasta), "application/gzip")
        return Response(archive_bytes.getvalue(), "application/gzip")

    monkeypatch.setattr(sources, "_open", open_source)
    releases = {
        "bold": sources.Release("BOLD_Public.09-Oct-2026", "https://doi.org/10.5883/DP-Latest", sources.BOLD_FASTA),
        "unite": sources.Release("2025-02-19", "https://doi.org/10.15156/BIO/3301232", "https://s3.hpc.ut.ee/unite.tgz"),
    }
    paths, provenance = sources.download_sources(releases, tmp_path)
    assert paths["bold"].read_bytes() == bold_fasta
    assert paths["unite"].read_bytes() == unite_fasta
    assert len(provenance["bold"]["sha256"]) == 64
    assert len(provenance["unite"]["sha256"]) == 64
    assert "uid" not in json.dumps(provenance)
    assert not (tmp_path / "unite.tgz").exists()


def test_ambiguous_bold_release_fails(monkeypatch):
    monkeypatch.setattr(sources, "_open", lambda *args, **kwargs: Response(b"<html></html>"))
    with pytest.raises(ValueError, match="BOLD Latest"):
        sources.discover_latest()
