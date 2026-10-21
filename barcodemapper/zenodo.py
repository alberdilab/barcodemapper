"""Resolve and download public Zenodo reference database records."""

import hashlib
import json
import os
import re
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import (DataHandler, FTPHandler, FileHandler,
                            HTTPRedirectHandler, Request, build_opener)

from barcodemapper.__version__ import __version__


ZENODO_URL = "https://zenodo.org"
DEFAULT_RECORD = None
USER_AGENT = f"BarcodeMapper/{__version__} (https://github.com/alberdilab/barcodemapper)"
METADATA_CAP = 1024 * 1024
DATABASE_CAP = 64 * 1024 ** 3


@dataclass(frozen=True)
class ZenodoFile:
    name: str
    size: int
    md5: str
    url: str
    _host: str | None = field(default=None, repr=False, compare=False)


@dataclass(frozen=True)
class ZenodoRecord:
    record_id: int
    doi: str
    version: str
    url: str
    license: str
    files: list[ZenodoFile]


class _SameHostRedirect(HTTPRedirectHandler):
    def __init__(self, host: str):
        self.host = host

    def redirect_request(self, request, response, code, message, headers, new_url):
        if _https_host(new_url) != self.host:
            raise ValueError("Zenodo redirect must use HTTPS on the base host")
        return super().redirect_request(request, response, code, message, headers, new_url)


class _NoFile(FileHandler):
    def file_open(self, request):
        raise ValueError("Non-HTTP Zenodo URL is forbidden")


class _NoFTP(FTPHandler):
    def ftp_open(self, request):
        raise ValueError("Non-HTTP Zenodo URL is forbidden")


class _NoData(DataHandler):
    def data_open(self, request):
        raise ValueError("Non-HTTP Zenodo URL is forbidden")


def _open(url: str, *, timeout: int = 60):
    host = _https_host(url)
    opener = build_opener(_SameHostRedirect(host), _NoFile(), _NoFTP(), _NoData())
    return opener.open(Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout)


def _https_host(url: str) -> str:
    parts = urlsplit(url)
    if parts.scheme != "https" or not parts.hostname or parts.username or parts.password:
        raise ValueError(f"Zenodo URL must use HTTPS with a valid host: {url!r}")
    return parts.hostname


def _base_url(base: str | None) -> str:
    url = (base if base is not None else os.environ.get("BARCODEMAPPER_ZENODO_URL") or ZENODO_URL).rstrip("/")
    _https_host(url)
    parts = urlsplit(url)
    if parts.query or parts.fragment:
        raise ValueError("Zenodo base URL cannot have a query or fragment")
    if parts.hostname not in ("zenodo.org", "sandbox.zenodo.org"):
        print(f"Warning: using nonstandard Zenodo host {parts.hostname}", file=sys.stderr)
    return url


def _file_name(name: str) -> str:
    if not isinstance(name, str) or not name or name.startswith(".") or ".." in name or any(char in name for char in ("/", "\\", "\x00")):
        raise ValueError(f"Unsafe Zenodo file name: {name!r}")
    return name


def _parse_record(data: dict, host: str) -> ZenodoRecord:
    if not isinstance(data.get("files"), list):
        raise ValueError("Zenodo record files must be a list")
    files = []
    for item in data["files"]:
        if not isinstance(item, dict) or not isinstance(item.get("links"), dict):
            raise ValueError("Zenodo record has invalid file metadata")
        name = _file_name(item.get("key", item.get("filename")))
        size = item.get("size", item.get("filesize"))
        if type(size) is not int or size < 0 or size > _size_cap(name):
            raise ValueError(f"Invalid Zenodo file size for {name!r}")
        url = item["links"].get("self", item["links"].get("download"))
        if _https_host(url) != host:
            raise ValueError(f"Zenodo download URL must use the base host: {url!r}")
        checksum = item["checksum"]
        if not isinstance(checksum, str) or not re.fullmatch(r"md5:[0-9a-fA-F]{32}", checksum):
            raise ValueError(f"Invalid Zenodo MD5 checksum for {name!r}")
        files.append(ZenodoFile(name, size, checksum[4:].lower(), url, host))
    license_id = data["metadata"]["license"]
    if isinstance(license_id, dict):
        license_id = license_id["id"]
    return ZenodoRecord(
        int(data["id"]), data["doi"], str(data["metadata"]["version"]),
        data["links"]["html"], license_id, files,
    )


def _size_cap(name: str) -> int:
    return DATABASE_CAP if name.endswith((".fa.gz", ".index.tar", ".tsv")) else METADATA_CAP


def resolve_record(record=None, version="latest", base=None) -> ZenodoRecord:
    """Resolve a concept or version record, optionally selecting a version string."""
    record = record if record is not None else DEFAULT_RECORD
    if record is None:
        raise ValueError("Zenodo record ID is required; DEFAULT_RECORD is not set")
    if not re.fullmatch(r"[0-9]+", str(record)):
        raise ValueError(f"Invalid Zenodo record ID: {record!r}")
    base_url = _base_url(base)
    host = _https_host(base_url)
    if version != "latest":
        page = f"{base_url}/api/records/{record}/versions?size=100"
        seen = set()
        while page:
            if page in seen or _https_host(page) != host:
                raise ValueError("Zenodo versions link is unsafe or repeated")
            seen.add(page)
            with _open(page) as response:
                data = json.load(response)
            match = next((item for item in data["hits"]["hits"]
                          if str(item.get("metadata", {}).get("version")) == version), None)
            if match is not None:
                record = match["id"]
                break
            page = data.get("links", {}).get("next")
        else:
            raise ValueError(f"Zenodo version {version!r} was not found for record {record}")
    with _open(f"{base_url}/api/records/{record}") as response:
        return _parse_record(json.load(response), host)


def download_file(file: ZenodoFile, directory, expected_sha256=None, progress=None) -> Path:
    """Download a file and report progress as (bytes_downloaded, total_bytes)."""
    name = _file_name(file.name)
    if type(file.size) is not int or file.size < 0 or file.size > _size_cap(name):
        raise ValueError(f"Invalid Zenodo file size for {name!r}")
    host = file._host or _https_host(_base_url(None))
    if _https_host(file.url) != host:
        raise ValueError(f"Zenodo download URL must use the base host: {file.url!r}")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    destination = directory / name
    temporary = None
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    downloaded = 0
    try:
        with _open(file.url, timeout=180) as response:
            with tempfile.NamedTemporaryFile(mode="wb", prefix=f".{name}.", suffix=".part", dir=directory, delete=False) as target:
                temporary = Path(target.name)
                for block in iter(lambda: response.read(1024 * 1024), b""):
                    target.write(block)
                    md5.update(block)
                    sha256.update(block)
                    downloaded += len(block)
                    if downloaded > file.size:
                        raise ValueError(f"Zenodo file exceeds declared size: {name}")
                    if progress is not None:
                        progress(downloaded, file.size)
        if downloaded != file.size:
            raise ValueError(f"Zenodo file size mismatch for {name}")
        if md5.hexdigest() != file.md5.lower():
            raise ValueError(f"Zenodo MD5 checksum mismatch for {name}")
        if expected_sha256 is not None and sha256.hexdigest() != expected_sha256.lower():
            raise ValueError(f"Zenodo SHA-256 checksum mismatch for {name}")
        temporary.replace(destination)
        return destination
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
