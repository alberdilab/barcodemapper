"""Check that a proposed release tag matches package metadata and changelog."""

import argparse
import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", help="Proposed release tag, for example v1.0.1")
    args = parser.parse_args()

    with (ROOT / "pyproject.toml").open("rb") as handle:
        version = tomllib.load(handle)["project"]["version"]

    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        parser.error(f"project version must be MAJOR.MINOR.PATCH: {version}")

    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    if not re.search(rf"^## \[{re.escape(version)}\](?:\s|$)", changelog, re.MULTILINE):
        parser.error(f"CHANGELOG.md needs a [{version}] section")

    if args.tag and args.tag != f"v{version}":
        parser.error(f"tag {args.tag} does not match project version v{version}")

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if not re.search(rf"^version:\s*[\"']?{re.escape(version)}[\"']?$", citation, re.MULTILINE):
        parser.error(f"CITATION.cff version must match {version}")

    if args.tag and not re.search(
        rf"^## \[{re.escape(version)}\] - \d{{4}}-\d{{2}}-\d{{2}}$",
        changelog,
        re.MULTILINE,
    ):
        parser.error(f"CHANGELOG.md needs a dated [{version}] release section before tagging")

    print(f"Version {version} matches CHANGELOG.md" + (f" and {args.tag}" if args.tag else ""))


if __name__ == "__main__":
    main()
