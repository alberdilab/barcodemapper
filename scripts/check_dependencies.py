"""Check pinned conda dependencies for newer releases."""

import argparse
import json
import re
import urllib.error
import urllib.request
from functools import cmp_to_key
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = re.compile(r"^(\s*-\s*)([A-Za-z0-9_.-]+)(\s*=\s*)([^\s#]+)")
NUMERIC_VERSION = re.compile(r"\d+(?:\.\d+)*\Z")


def parse_environment(text: str) -> tuple[list[str], dict[str, str]]:
    channels = []
    pins = {}
    section = None

    for line in text.splitlines():
        if line and not line[0].isspace() and ":" in line:
            section, value = line.split(":", 1)
            if section == "channels" and value.strip().startswith("["):
                channels.extend(item.strip().strip("'\"") for item in value.strip()[1:-1].split(","))
            continue
        if section == "channels":
            match = re.match(r"^\s*-\s*([^\s#]+)", line)
            if match:
                channels.append(match.group(1).strip("'\""))
        elif section == "dependencies":
            match = DEPENDENCY.match(line)
            if match:
                pins[match.group(2)] = match.group(4)

    return channels, pins


def compare_versions(left: str, right: str) -> int:
    if NUMERIC_VERSION.fullmatch(left) and NUMERIC_VERSION.fullmatch(right):
        left_parts = tuple(int(part) for part in left.split("."))
        right_parts = tuple(int(part) for part in right.split("."))
        length = max(len(left_parts), len(right_parts))
        left_parts += (0,) * (length - len(left_parts))
        right_parts += (0,) * (length - len(right_parts))
        return (left_parts > right_parts) - (left_parts < right_parts)
    return (left.lower() > right.lower()) - (left.lower() < right.lower())


def fetch_latest(channel: str, name: str) -> str | None:
    url = f"https://api.anaconda.org/package/{channel}/{name}"
    request = urllib.request.Request(url, headers={"User-Agent": "barcodemapper-dependency-check/1.0"})
    for attempt in range(2):
        try:
            with urllib.request.urlopen(request, timeout=10) as response:
                version = json.load(response)["latest_version"]
            if not isinstance(version, str) or not version:
                raise ValueError(f"invalid latest_version for {channel}/{name}")
            return version
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return None
            if attempt:
                raise
        except (OSError, ValueError, KeyError, TypeError):
            if attempt:
                raise
    return None


def write_candidate(path: Path, source: str, versions: dict[str, str], bump_python: bool, exact: bool = False) -> None:
    """Write the candidate env file.

    By default every pin except python's minor is dropped so the solver picks the
    newest versions that work together; with exact=True each pin becomes the
    latest release of that package (which may not be installable as a set).
    """
    lines = source.splitlines(keepends=True)
    section = None
    for index, line in enumerate(lines):
        if line.startswith("name:"):
            ending = "\n" if line.endswith("\n") else ""
            lines[index] = "name: barcodemapper-candidate" + ending
        if line and not line[0].isspace() and ":" in line:
            section = line.split(":", 1)[0]
            continue
        if section != "dependencies":
            continue
        match = DEPENDENCY.match(line)
        if not match or match.group(2) not in versions:
            continue
        if match.group(2) == "python" and not bump_python:
            continue
        if exact:
            lines[index] = line[:match.start(4)] + versions[match.group(2)] + line[match.end(4):]
        else:
            lines[index] = line[:match.end(2)] + line[match.end(4):]
    path.write_text("".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-candidate", type=Path, metavar="PATH", help="Write an updated candidate env file")
    parser.add_argument("--exact", action="store_true",
                        help="Pin the candidate to each package's latest release instead of leaving the solver to "
                             "pick the newest compatible set")
    parser.add_argument("--bump-python", action="store_true", help="Update the Python pin in the candidate")
    parser.add_argument("--json", action="store_true", help="Print machine-readable results")
    parser.add_argument("--fail-if-outdated", action="store_true", help="Exit 10 if any dependency is outdated")
    args = parser.parse_args(argv)

    source = (ROOT / "barcodemapper_env.yml").read_text(encoding="utf-8")
    channels, pins = parse_environment(source)
    results = {"outdated": [], "current": [], "failed": []}
    rows = []
    versions = {}
    for name, pinned in pins.items():
        try:
            available = [version for channel in channels if (version := fetch_latest(channel, name)) is not None]
            if not available:
                raise ValueError(f"no releases found for {name}")
            latest = max(available, key=cmp_to_key(compare_versions))
            category = "outdated" if compare_versions(latest, pinned) > 0 else "current"
            held = name == "python" and not args.bump_python
            if held:
                category = "current"  # the Python minor stays pinned unless --bump-python
            versions[name] = latest
        except (OSError, ValueError, KeyError, TypeError):
            latest = None
            category = "failed"
            held = False
        row = {"package": name, "pinned": pinned, "latest": latest,
               "status": "held (--bump-python)" if held and compare_versions(latest, pinned) > 0 else
                    {"outdated": "newer available", "current": "current", "failed": "lookup failed"}[category]}
        results[category].append(row)
        rows.append(row)

    if args.write_candidate:
        write_candidate(args.write_candidate, source, versions, args.bump_python, args.exact)
    if args.json:
        print(json.dumps(results))
    else:
        widths = [max(len(str(row[key] or "-")) for row in rows) if rows else 0
                  for key in ("package", "pinned", "latest")]
        widths = [max(width, len(label)) for width, label in zip(widths, ("package", "pinned", "latest"))]
        print(f"{'package':<{widths[0]}}  {'pinned':<{widths[1]}}  {'latest':<{widths[2]}}  status")
        for row in rows:
            print(f"{row['package']:<{widths[0]}}  {row['pinned']:<{widths[1]}}  "
                  f"{(row['latest'] or '-'):<{widths[2]}}  {row['status']}")

    if results["failed"]:
        return 2
    if args.fail_if_outdated and results["outdated"]:
        return 10
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
