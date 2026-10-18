"""Run the small end-to-end fixture and compare its taxonomy table to the golden."""

import argparse
import difflib
import importlib.metadata
import json
import sys
import tempfile
from contextlib import nullcontext
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "regression" / "data"
GOLDEN = ROOT / "tests" / "regression" / "expected_taxonomy.tsv"


def _sorted_rows(table: str) -> tuple[str, int, int]:
    lines = table.splitlines(keepends=True)
    if not lines or not lines[0].startswith("Taxonomy\t"):
        raise ValueError("taxonomy table has no Taxonomy header")
    return lines[0] + "".join(sorted(lines[1:])), len(lines) - 1, len(lines[0].split("\t")) - 1


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _print_versions(manifest: dict) -> None:
    for tool, version in sorted(manifest["tools"].items()):
        print(f"{tool}: {version}")
    print(f"python: {sys.version.split()[0]}")
    print(f"pandas: {_version('pandas')}")
    print(f"pysam: {_version('pysam')}")


def _run(executable: str, directory: Path, update_golden: bool) -> int:
    output = directory / "taxonomy.tsv"
    command = [
        executable, "-i", str(DATA / "reads"), "-d", str(DATA / "references.fa"),
        "-o", str(output), "-t", str(directory / "work"),
    ]
    try:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    except OSError as exc:
        print(f"regression check: pipeline failed: {exc}", file=sys.stderr)
        return 2
    if result.returncode != 0:
        print(f"regression check: pipeline failed (exit {result.returncode})", file=sys.stderr)
        print("\n".join(result.stderr.splitlines()[-40:]), file=sys.stderr)
        return 2

    manifest_path = output.with_suffix(".manifest.json")
    try:
        actual, taxa, samples = _sorted_rows(output.read_bytes().decode("utf-8"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        _print_versions(manifest)
    except (OSError, UnicodeError, ValueError, KeyError, AttributeError) as exc:
        print(f"regression check: pipeline output is incomplete: {exc}", file=sys.stderr)
        return 2

    if taxa == 0 or samples == 0:
        print("regression check: the pipeline produced an empty taxonomy table", file=sys.stderr)
        return 2
    if update_golden:
        GOLDEN.write_text(actual, encoding="utf-8", newline="")
        print(f"regression check: updated {GOLDEN} ({taxa} taxa, {samples} samples)")
        return 0

    try:
        expected, _, _ = _sorted_rows(GOLDEN.read_bytes().decode("utf-8"))
    except FileNotFoundError:
        expected = ""
        print(f"regression check: golden file missing: {GOLDEN}", file=sys.stderr)
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"regression check: invalid golden file: {exc}", file=sys.stderr)
        return 1
    if actual != expected:
        sys.stdout.writelines(difflib.unified_diff(
            expected.splitlines(keepends=True), actual.splitlines(keepends=True),
            fromfile=str(GOLDEN), tofile=str(output),
        ))
        return 1
    print(f"regression check: OK ({taxa} taxa, {samples} samples)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--barcodemapper", default="barcodemapper", help="Installed barcodemapper entry point")
    parser.add_argument("--keep-dir", action="store_true", help="Keep the temporary run directory")
    parser.add_argument("--update-golden", action="store_true", help="Replace the committed taxonomy table")
    args = parser.parse_args()

    directory = (nullcontext(Path(tempfile.mkdtemp(prefix="bm-regression-"))) if args.keep_dir
                 else tempfile.TemporaryDirectory(prefix="bm-regression-"))
    with directory as path:
        path = Path(path)
        if args.keep_dir:
            print(f"regression check: run directory {path}")
        return _run(args.barcodemapper, path, args.update_golden)


if __name__ == "__main__":
    raise SystemExit(main())
