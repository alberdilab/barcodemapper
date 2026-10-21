"""Prepare portable database manifests and publication credit."""

from copy import deepcopy
from getpass import getuser
from html import escape
from pathlib import Path, PurePosixPath, PureWindowsPath
import re


BOLD_V4 = ("Ratnasingham S, Wei C, Chan D, et al. (2024) BOLD v4: A Centralized "
           "Bioinformatics Platform for DNA-Based Biodiversity Data. In: DNA Barcoding: "
           "Methods and Protocols, pp. 403-441. Springer US")
BOLD_ALIAS = "10.5883/DP-Latest is a weekly-updating alias (moving alias)"
CC_BY_SA = "https://creativecommons.org/licenses/by-sa/4.0/"
CC_BY = "https://creativecommons.org/licenses/by/4.0/"
_PATH = re.compile(r"(?:/projects|/home|/Users)/[^\s,;)'\"\]}]+|[A-Za-z]:\\[^\s,;)'\"\]}]+|~[/\\][^\s,;)'\"\]}]+")


def _is_absolute_path(value: str) -> bool:
    return PurePosixPath(value).is_absolute() or PureWindowsPath(value).is_absolute()


def _scrub_paths(value):
    if isinstance(value, dict):
        return {key: _scrub_paths(item) for key, item in value.items()
                if not (isinstance(key, str) and (_is_absolute_path(key) or _PATH.search(key)))}
    if isinstance(value, list):
        return [_scrub_paths(item) for item in value]
    if isinstance(value, str):
        if _is_absolute_path(value):
            return ""
        value = value.replace(str(Path.home()), "[redacted]")
        value = _PATH.sub("[redacted]", value)
        user = getuser()
        if user:
            value = re.sub(rf"(?<![\w]){re.escape(user)}(?![\w])", "[redacted]", value)
    return value


def validate_source_manifest(manifest: dict, fasta_hash: str) -> None:
    """Require an automatically built database with usable source provenance."""
    if not isinstance(manifest, dict) or not isinstance(manifest.get("database"), dict):
        raise ValueError("database manifest is invalid")
    if manifest["database"].get("sha256") != fasta_hash:
        raise ValueError("manifest FASTA SHA-256 does not match the database")
    sources, parameters = manifest.get("sources"), manifest.get("parameters")
    if manifest.get("automatic_sources") is not True or not isinstance(sources, dict) or not isinstance(parameters, dict):
        raise ValueError("manifest lacks automatic source provenance")
    for name in ("bold", "unite"):
        source = sources.get(name)
        if not isinstance(source, dict) or not all(source.get(key) for key in ("release", "doi")) or not parameters.get(f"{name}_retain"):
            raise ValueError(f"manifest lacks {name.upper()} release, DOI or retained-taxa provenance")


def portable_manifest(
    manifest: dict, *, origin: str = "zenodo", record_doi: str | None = None,
    record_url: str | None = None,
) -> dict:
    """Copy a manifest while removing local absolute paths."""
    portable = _scrub_paths(deepcopy(manifest))
    parameters = manifest.get("parameters")
    if isinstance(parameters, dict):
        portable["parameters"] = {
            key: _scrub_paths(value) for key, value in parameters.items()
            if isinstance(key, str) and not re.search(r"path|dir|db|file|output|tmp", key, re.I)
            and isinstance(value, (str, int, float, bool, type(None)))
        }
    tools = manifest.get("tools")
    if isinstance(tools, dict):
        portable["tools"] = {key: _scrub_paths(value) for key, value in tools.items()
                             if isinstance(key, str) and not _PATH.search(key) and isinstance(value, str)}
    database = portable.get("database", {})
    original_path = manifest.get("database", {}).get("path")
    if isinstance(original_path, str):
        database["path"] = PureWindowsPath(original_path).name
    portable["origin"] = origin
    if record_doi is not None or record_url is not None:
        portable["zenodo"] = {"doi": record_doi, "url": record_url}
    return portable


def _source_description(name: str, source: dict) -> str:
    release = source.get("release")
    doi = source.get("doi")
    digest = source.get("sha256")
    if release:
        label = f"{name} release {release}"
    else:
        label = f"{name} user-supplied source"
    details = []
    if doi:
        details.append(f"DOI: {doi}")
    if digest:
        details.append(f"SHA-256: {digest}")
    return escape(label + (f" ({'; '.join(details)})" if details else ""))


def _unite_citation(source: dict) -> str:
    release = str(source.get("release") or "")
    year = re.search(r"(?:19|20)\d{2}", release)
    year_text = f" ({year.group()})" if year else ""
    return ("Abarenkov K, Zirk A, Piirmann T, Pöhönen R, Ivanov F, Nilsson RH, "
            f"Kõljalg U{year_text}: UNITE general FASTA release for eukaryotes 2. "
            f"UNITE Community. {source.get('doi') or 'not recorded'}")


def zenodo_metadata(
    manifest: dict, *, version: str, creators: list[dict],
    github_url: str = "https://github.com/alberdilab/barcodemapper",
) -> dict:
    """Return Zenodo deposit metadata for a combined reference database."""
    sources = manifest.get("sources") or {}
    parameters = manifest.get("parameters") or {}
    bold = sources.get("bold") or {}
    unite = sources.get("unite") or {}
    description = (
        f"<p>Combined BOLD and UNITE reference database derived from "
        f"{_source_description('BOLD', bold)} and {_source_description('UNITE', unite)}.</p>"
        f"<p>The source data were modified: merged, dereplicated, filtered to the retained taxa "
        f"(BOLD: {escape(str(parameters.get('bold_retain') or 'not recorded'))}; "
        f"UNITE: {escape(str(parameters.get('unite_retain') or 'not recorded'))}), "
        f"and sequence headers reformatted. The database is released under CC BY-SA 4.0. "
        f"BarcodeMapper code is licensed under MIT.</p>"
        f"<p>Accessed/built: {escape(str(manifest.get('created_utc') or '')[:10])}. "
        f"BOLD snapshot release: {escape(str(bold.get('release') or 'not recorded'))}. "
        f"{escape(BOLD_ALIAS)}. BOLD source: CC BY-SA 4.0 ({CC_BY_SA}). "
        f"The UNITE-derived portion is also available from UNITE under CC BY 4.0 ({CC_BY}). "
        f"Whole database: CC BY-SA 4.0 ({CC_BY_SA}).</p>"
        f"<p>Please cite BOLD: Ratnasingham &amp; Hebert 2007 Mol Ecol Notes 7:355-364; "
        f"Ratnasingham &amp; Hebert 2013 PLoS ONE 8(8):e66213; {escape(BOLD_V4)}.</p>"
        f"<p>Please cite UNITE: {escape(_unite_citation(unite))}.</p>"
    )
    related = []
    for source in (bold, unite):
        doi = source.get("doi")
        if doi:
            related.append({"identifier": doi.removeprefix("https://doi.org/"),
                            "relation": "isDerivedFrom", "resource_type": "dataset", "scheme": "doi"})
    related.append({"identifier": github_url, "relation": "isSupplementTo",
                    "resource_type": "software", "scheme": "url"})
    return {
        "title": "BarcodeMapper combined BOLD and UNITE reference database",
        "upload_type": "dataset",
        "version": version,
        "license": "cc-by-sa-4.0",
        "creators": creators,
        "description": description,
        "related_identifiers": related,
    }


def attribution_text(manifest: dict) -> str:
    """Return plain-text licence, source provenance, and citation guidance."""
    sources = manifest.get("sources") or {}
    parameters = manifest.get("parameters") or {}
    lines = [
        "BarcodeMapper combined BOLD and UNITE reference database",
        "",
        f"Licence: CC BY-SA 4.0 ({CC_BY_SA}).",
        f"Accessed/built: {str(manifest.get('created_utc') or '')[:10]}",
        f"BOLD snapshot release: {(sources.get('bold') or {}).get('release') or 'not recorded'}",
        f"{BOLD_ALIAS}.",
        f"BOLD: CC BY-SA 4.0 ({CC_BY_SA}).",
        f"The UNITE-derived portion is also available from UNITE under CC BY 4.0 ({CC_BY}).",
        "",
    ]
    for name in ("bold", "unite"):
        source = sources.get(name) or {}
        lines.extend([
            f"{name.upper()} source:",
            f"  Release: {source.get('release') or 'user-supplied'}",
            f"  DOI: {source.get('doi') or 'not recorded'}",
            f"  Source URL: {source.get('url') or 'not recorded'}",
            f"  SHA-256: {source.get('sha256') or 'not recorded'}",
            "",
        ])
    lines.extend([
        "Changes made: The BOLD and UNITE sources were merged, dereplicated, filtered to",
        f"retained taxa (BOLD: {parameters.get('bold_retain') or 'not recorded'}; "
        f"UNITE: {parameters.get('unite_retain') or 'not recorded'}), and sequence headers",
        "were reformatted.",
        "",
        "Please cite BOLD: Ratnasingham & Hebert 2007 Mol Ecol Notes 7:355-364;",
        "Ratnasingham & Hebert 2013 PLoS ONE 8(8):e66213.",
        BOLD_V4,
        f"Please cite UNITE: {_unite_citation(sources.get('unite') or {})}.",
    ])
    return "\n".join(lines) + "\n"
