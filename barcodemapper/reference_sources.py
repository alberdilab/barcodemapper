"""Discover and download the current BOLD and UNITE reference releases."""

import gzip
import hashlib
import json
import re
import shutil
import tarfile
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen

from barcodemapper.__version__ import __version__


BOLD_LATEST = "https://bench.boldsystems.org/index.php/datapackages/Latest"
BOLD_FASTA = "https://bench.boldsystems.org/index.php/API_Datapackage/fasta"
UNITE_REPOSITORY = "https://unite.ut.ee/repository.php"
PLUTOF_DOIS = "https://api.plutof.ut.ee/v1/public/dois/"
USER_AGENT = f"BarcodeMapper/{__version__} (https://github.com/alberdilab/barcodemapper)"


@dataclass(frozen=True)
class Release:
    identifier: str
    doi: str
    url: str


def _open(url: str, *, timeout: int = 60):
    return urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout)


class _BoldPackageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.identifiers = []

    def handle_starttag(self, tag, attrs):
        if tag != "button":
            return
        attrs = dict(attrs)
        if attrs.get("name") in {"download-fasta", "download-fasta-unauthenticated"}:
            self.identifiers.append(attrs.get("data-package-id"))


class _UniteReleasesParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.next_table = False
        self.in_table = False
        self.row = None
        self.cell = None
        self.rows = []

    def handle_starttag(self, tag, attrs):
        if tag == "table" and self.next_table:
            self.in_table = True
            self.next_table = False
        elif self.in_table and tag == "tr":
            self.row = []
        elif self.in_table and tag == "td":
            self.cell = {"text": "", "href": None}
        elif self.cell is not None and tag == "a":
            self.cell["href"] = dict(attrs).get("href")

    def handle_data(self, data):
        if "List of all general FASTA releases" in data:
            self.next_table = True
        if self.cell is not None:
            self.cell["text"] += data

    def handle_endtag(self, tag):
        if tag == "td" and self.cell is not None:
            self.cell["text"] = " ".join(self.cell["text"].split())
            self.row.append(self.cell)
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.rows.append(self.row)
            self.row = None
        elif tag == "table" and self.in_table:
            self.in_table = False


def discover_latest() -> dict[str, Release]:
    """Resolve the current BOLD snapshot and UNITE general FASTA with singletons."""
    with _open(BOLD_LATEST) as response:
        bold_page = response.read().decode("utf-8")
    bold_parser = _BoldPackageParser()
    bold_parser.feed(bold_page)
    bold_ids = set(bold_parser.identifiers)
    if len(bold_ids) != 1:
        raise ValueError("Could not identify one BOLD Latest FASTA package")
    bold_id = bold_ids.pop()
    if not bold_id or not re.fullmatch(r"BOLD_Public\.\d{2}-[A-Za-z]{3}-\d{4}", bold_id):
        raise ValueError(f"Unexpected BOLD package identifier: {bold_id!r}")

    with _open(UNITE_REPOSITORY) as response:
        unite_page = response.read().decode("utf-8")
    unite_parser = _UniteReleasesParser()
    unite_parser.feed(unite_page)
    candidates = []
    for row in unite_parser.rows:
        if len(row) < 8 or row[2]["text"].casefold() != "all eukaryotes" or row[5]["text"].casefold() != "current":
            continue
        try:
            date = row[1]["text"]
            count = int(row[4]["text"].replace(" ", ""))
            doi = row[6]["href"]
            if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date) or not doi or not re.fullmatch(r"https://doi\.org/10\.15156/BIO/\d+", doi):
                continue
            candidates.append((date, count, doi))
        except (ValueError, TypeError):
            continue
    if not candidates:
        raise ValueError("Could not identify the current UNITE All Eukaryotes general FASTA release")
    # The larger current All Eukaryotes release includes the singleton sequences
    # used by BarcodeMapper's earlier database snapshot.
    date, _, doi = max(candidates)
    identifier = doi.removeprefix("https://doi.org/")
    query = urlencode({"format": "vnd.api+json", "identifier": identifier})
    with _open(f"{PLUTOF_DOIS}?{query}") as response:
        data = json.load(response)
    entries = data.get("data", [])
    if len(entries) != 1:
        raise ValueError(f"PlutoF did not return one record for {identifier}")
    media = entries[0].get("attributes", {}).get("media", [])
    archives = [item for item in media if re.search(r"general_release_(?:dynamic_)?s_all_.*\.tgz$", item.get("name", ""))]
    if not archives:
        raise ValueError(f"PlutoF did not provide a matching UNITE archive for {identifier}")
    unite_url = archives[-1].get("url", "")
    if urlparse(unite_url).scheme != "https":
        raise ValueError("UNITE archive URL must use HTTPS")
    return {
        "bold": Release(bold_id, "https://doi.org/10.5883/DP-Latest", BOLD_FASTA),
        "unite": Release(date, doi, unite_url),
    }


def _copy_and_hash(source, destination: Path) -> str:
    digest = hashlib.sha256()
    with destination.open("wb") as target:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            target.write(block)
            digest.update(block)
    if destination.stat().st_size == 0:
        raise ValueError(f"Downloaded reference is empty: {destination.name}")
    return digest.hexdigest()


def download_sources(releases: dict[str, Release], directory: Path) -> tuple[dict[str, Path], dict[str, dict]]:
    """Download and extract the two FASTAs, returning paths and provenance."""
    directory.mkdir(parents=True, exist_ok=True)
    bold = releases["bold"]
    bold_query = urlencode({"id": bold.identifier})
    with _open(f"{bold.url}?{bold_query}") as response:
        uid = json.load(response)
    if not isinstance(uid, str) or not re.fullmatch(r"[A-Za-z0-9]+", uid):
        raise ValueError("BOLD did not generate a valid FASTA download link")
    bold_path = directory / "bold.fasta"
    with _open(f"{bold.url}?{urlencode({'id': bold.identifier, 'uid': uid})}", timeout=180) as response:
        if response.headers.get_content_type() not in {"application/x-gzip", "application/gzip", "application/octet-stream"}:
            raise ValueError("BOLD returned a response other than a compressed FASTA")
        with gzip.GzipFile(fileobj=response) as uncompressed:
            bold_sha = _copy_and_hash(uncompressed, bold_path)

    unite = releases["unite"]
    archive_path = directory / "unite.tgz"
    with _open(unite.url, timeout=180) as response:
        _copy_and_hash(response, archive_path)
    unite_path = directory / "unite.fasta"
    with tarfile.open(archive_path, "r:gz") as archive:
        members = [member for member in archive.getmembers() if
                   member.isfile() and not member.name.endswith("_dev.fasta") and
                   re.fullmatch(r"(?:.*/)?sh_general_release_dynamic_s_all_.*\.fasta", member.name)]
        if len(members) != 1:
            raise ValueError("UNITE archive must contain one standard dynamic All Eukaryotes FASTA")
        with archive.extractfile(members[0]) as source:
            unite_sha = _copy_and_hash(source, unite_path)
    archive_path.unlink()
    paths = {"bold": bold_path, "unite": unite_path}
    provenance = {
        "bold": {"release": bold.identifier, "doi": bold.doi, "url": f"{bold.url}?{bold_query}", "sha256": bold_sha},
        "unite": {"release": unite.identifier, "doi": unite.doi, "url": unite.url, "sha256": unite_sha},
    }
    return paths, provenance
