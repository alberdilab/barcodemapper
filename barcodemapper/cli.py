"""Command line interface for the BarcodeMapper workflow."""

import argparse
import hashlib
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import yaml

from barcodemapper.__version__ import __version__
from barcodemapper.reference_sources import Release, discover_latest, download_sources
from barcodemapper.workflow.scripts.validate_database import validate_database


PACKAGE_DIR = Path(__file__).resolve().parent
READ_PATTERN = re.compile(r"^(?P<sample>.+)_(?P<mate>[12])\.(?:fq|fastq)\.gz$")
DEFAULT_ADAPTER = "CTGTCTCTTATACACATCT"
DEFAULT_MAX_HITS = 200
RESERVED_SAMPLE_NAMES = {"Taxonomy", "Phylum", "Kingdom", "color"}
PRESETS = {
    "ultrastrict": (0.0, 120),
    "strict": (2.0, 100),
    "standard": (3.0, 80),
    "permissive": (5.0, 60),
}
INDEX_PARTS = ("1", "2", "3", "4", "rev.1", "rev.2")


def _database_base(path: Path) -> Path:
    """Return the artifact stem for plain and gzip-compressed FASTA paths."""
    return path.with_suffix("").with_suffix("") if path.suffix.lower() == ".gz" else path.with_suffix("")


def _database_artifact(path: Path, suffix: str) -> Path:
    return Path(f"{_database_base(path)}{suffix}")


def _current_files(paths: list[Path], database: Path) -> bool:
    database_mtime = database.stat().st_mtime_ns
    return all(
        path.is_file() and path.stat().st_size > 0 and path.stat().st_mtime_ns >= database_mtime
        for path in paths
    )


def _database_artifacts(database: Path, work_dir: Path, build: bool) -> tuple[Path, str, bool, Path, bool]:
    """Choose existing artifacts or place generated analysis artifacts in the work directory."""
    database_index = _database_base(database)
    database_headers = _database_artifact(database, ".tsv")
    if build:
        return database_index, "bt2l", True, database_headers, True

    for extension in ("bt2l", "bt2"):
        index_files = [Path(f"{database_index}.{part}.{extension}") for part in INDEX_PARTS]
        if _current_files(index_files, database):
            index_base, generate_index = database_index, False
            break
    else:
        index_base = work_dir / "reference" / database_index.name
        extension, generate_index = "bt2l", True

    if _current_files([database_headers], database):
        headers_path, generate_headers = database_headers, False
    else:
        headers_path, generate_headers = work_dir / "reference" / f"{database_index.name}.tsv", True
    return index_base, extension, generate_index, headers_path, generate_headers


def _snakemake_command(tmp_dir: Path, config_path: Path, profile: str) -> list[str]:
    return [
        "snakemake",
        "-s", str(PACKAGE_DIR / "workflow" / "Snakefile"),
        "--directory", str(tmp_dir),
        "--workflow-profile", str(PACKAGE_DIR / "profile" / profile),
        "--configfile", str(config_path),
        "--config", f"package_dir={PACKAGE_DIR}",
    ]


def unlock_snakemake(tmp_dir: Path, config_path: Path, profile: str) -> None:
    subprocess.run(_snakemake_command(tmp_dir, config_path, profile) + ["--unlock"], check=True)
    print(f"The working directory {tmp_dir} has been unlocked.")


def run_snakemake(tmp_dir: Path, config_path: Path, profile: str) -> None:
    subprocess.run(_snakemake_command(tmp_dir, config_path, profile), check=True)


def _sample_maps(paths: list[Path]) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    mates: dict[int, dict[str, list[str]]] = {1: defaultdict(list), 2: defaultdict(list)}
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Read file does not exist: {path}")
        match = READ_PATTERN.fullmatch(path.name)
        if match is None:
            raise ValueError(f"Read filename must end in _1/_2.fq.gz or _1/_2.fastq.gz: {path.name}")
        sample = match.group("sample")
        if sample in RESERVED_SAMPLE_NAMES:
            raise ValueError(f"Sample name {sample!r} is reserved for an output column")
        mate = int(match.group("mate"))
        mates[mate][sample].append(str(path.resolve()))

    forward = mates[1]
    reverse = mates[2]
    if not forward and not reverse:
        raise ValueError("No paired FASTQ files were found")
    if set(forward) != set(reverse):
        missing_r1 = sorted(set(reverse) - set(forward))
        missing_r2 = sorted(set(forward) - set(reverse))
        raise ValueError(f"Unpaired samples: missing R1={missing_r1}, missing R2={missing_r2}")
    for sample in sorted(forward):
        if len(forward[sample]) != 1 or len(reverse[sample]) != 1:
            raise ValueError(
                f"Sample {sample!r} has multiple files per mate. Merge lanes first or give each pair a distinct sample name."
            )
    return dict(sorted(forward.items())), dict(sorted(reverse.items()))


def _save_samples(
    read1: dict[str, list[str]], read2: dict[str, list[str]], output_dir: Path
) -> None:
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    for mate, mapping in ((1, read1), (2, read2)):
        (data_dir / f"sample_to_reads{mate}.json").write_text(
            json.dumps(mapping, indent=2) + "\n", encoding="utf-8"
        )


def inputdir_to_samples(input_dir: str | Path, output_dir: str | Path) -> tuple[dict, dict]:
    directory = Path(input_dir).expanduser().resolve()
    if not directory.is_dir():
        raise ValueError(f"Input directory does not exist: {directory}")
    paths = sorted(
        path for path in directory.iterdir()
        if path.name.endswith((".fq.gz", ".fastq.gz"))
    )
    maps = _sample_maps(paths)
    _save_samples(*maps, Path(output_dir))
    return maps


def inputlist_to_samples(read1: str, read2: str, output_dir: str | Path) -> tuple[dict, dict]:
    def split_files(value: str, mate: int) -> list[Path]:
        paths = [Path(item.strip()).expanduser().resolve() for item in value.split(",") if item.strip()]
        if not paths:
            raise ValueError(f"No R{mate} files were supplied")
        for path in paths:
            match = READ_PATTERN.fullmatch(path.name)
            if match is None or int(match.group("mate")) != mate:
                raise ValueError(f"Expected an R{mate} filename: {path}")
        return paths

    maps = _sample_maps(split_files(read1, 1) + split_files(read2, 2))
    _save_samples(*maps, Path(output_dir))
    return maps


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tool_version(name: str) -> str:
    try:
        completed = subprocess.run([name, "--version"], capture_output=True, text=True, check=False)
    except OSError:
        return "unavailable"
    lines = (completed.stdout + completed.stderr).splitlines()
    return lines[0].strip() if lines else "unknown"


def _write_manifest(
    destination: Path,
    db_path: Path,
    config: dict,
    read1: dict[str, list[str]],
    read2: dict[str, list[str]],
    source_info: dict | None = None,
    database_path: Path | None = None,
    automatic_sources: bool = False,
) -> None:
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "barcodemapper_version": __version__,
        "database": {"path": str(database_path or db_path), "sha256": _sha256(db_path)},
        "parameters": config,
        "samples": {sample: {"read1": read1[sample], "read2": read2[sample]} for sample in read1},
        "tools": {tool: _tool_version(tool) for tool in ("snakemake", "fastp", "bowtie2", "bowtie2-build", "samtools")},
    }
    if automatic_sources:
        manifest["automatic_sources"] = True
    if source_info is not None:
        manifest["sources"] = source_info
    elif "bold_db" in config:
        manifest["sources"] = {
            source: {"path": config[key], "sha256": _sha256(Path(config[key]))}
            for source, key in (("bold", "bold_db"), ("unite", "unite_db"))
        }
    for sample in read1:
        fastp_file = Path(config["tmp_dir"]) / "fastp" / f"{sample}.json"
        taxonomy_file = Path(config["tmp_dir"]) / "taxonomy" / f"{sample}.tsv"
        if fastp_file.is_file():
            fastp_data = json.loads(fastp_file.read_text(encoding="utf-8"))
            manifest["samples"][sample]["fastp_reads_before"] = fastp_data.get("summary", {}).get("before_filtering", {}).get("total_reads")
            manifest["samples"][sample]["fastp_reads_after"] = fastp_data.get("summary", {}).get("after_filtering", {}).get("total_reads")
        if taxonomy_file.is_file():
            with taxonomy_file.open(encoding="utf-8") as handle:
                manifest["samples"][sample]["assigned_read_ids"] = max(0, sum(1 for _ in handle) - 1)
    destination.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _auto_manifest(db_path: Path) -> dict | None:
    path = _database_artifact(db_path, ".manifest.json")
    if not path.is_file():
        return None
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return manifest if isinstance(manifest, dict) and manifest.get("automatic_sources") is True else None


def _same_auto_database(db_path: Path, manifest: dict | None, releases: dict[str, Release], args) -> bool:
    if not db_path.is_file() or not manifest:
        return False
    params = manifest.get("parameters")
    sources = manifest.get("sources")
    database = manifest.get("database")
    if not all(isinstance(value, dict) for value in (params, sources, database)):
        return False
    if params.get("bold_retain") != args.bold_retain:
        return False
    if params.get("unite_retain") != args.unite_retain:
        return False
    for name in ("bold", "unite"):
        source = sources.get(name, {})
        if not isinstance(source, dict):
            return False
        if source.get("release") != releases[name].identifier or source.get("doi") != releases[name].doi:
            return False
    return database.get("sha256") == _sha256(db_path)


def _build_config(db_path: Path, tmp_dir: Path, args, *, bold: Path | None = None,
                  unite: Path | None = None, generate_index: bool = True,
                  generate_headers: bool = True) -> dict:
    config = {
        "run": datetime.now(timezone.utc).isoformat(),
        "version": __version__,
        "barcodemapper_db": str(db_path),
        "barcodemapper_index_base": str(_database_base(db_path)),
        "barcodemapper_index_extension": "bt2l",
        "generate_index": generate_index,
        "barcodemapper_db_headers": str(_database_artifact(db_path, ".tsv")),
        "generate_headers": generate_headers,
        "tmp_dir": str(tmp_dir),
        "output_file": "no",
        "build_only": "yes",
        "report_file": "none",
    }
    if bold is not None and unite is not None:
        config.update({"bold_db": str(bold), "unite_db": str(unite),
                       "bold_retain": args.bold_retain, "unite_retain": args.unite_retain})
    return config


def _promote_database(staged: Path, destination: Path) -> None:
    """Replace database artifacts after all new artifacts are complete."""
    base = _database_base(staged)
    names = [staged.name, _database_artifact(staged, ".tsv").name]
    names += [f"{base.name}.{part}.bt2l" for part in INDEX_PARTS]
    names.append(_database_artifact(staged, ".manifest.json").name)
    old_names = names + [f"{base.name}.{part}.bt2" for part in INDEX_PARTS]
    backup = staged.parent / "previous"
    if backup.exists():
        shutil.rmtree(backup)
    backup.mkdir()
    moved_old = []
    moved_new = []
    try:
        for name in old_names:
            old = destination.parent / name
            if old.exists():
                old.replace(backup / name)
                moved_old.append(name)
        for name in names:
            (staged.parent / name).replace(destination.parent / name)
            moved_new.append(name)
    except OSError as exc:
        rollback_failed = False
        try:
            for name in reversed(moved_new):
                (destination.parent / name).replace(staged.parent / name)
            for name in reversed(moved_old):
                (backup / name).replace(destination.parent / name)
        except OSError as rollback_exc:
            rollback_failed = True
            exc.add_note(f"database promotion rollback failed: {rollback_exc}")
        finally:
            if not rollback_failed:
                try:
                    shutil.rmtree(backup)
                except OSError as cleanup_exc:
                    exc.add_note(f"database promotion cleanup failed: {cleanup_exc}")
        raise
    shutil.rmtree(backup)


def _build_outputs(db_path: Path) -> list[Path]:
    base = _database_base(db_path)
    return [db_path, _database_artifact(db_path, ".tsv")] + [
        Path(f"{base}.{part}.bt2l") for part in INDEX_PARTS
    ]


def _manual_target_guard(parser: argparse.ArgumentParser, db_path: Path) -> None:
    manifest = _database_artifact(db_path, ".manifest.json")
    base = _database_base(db_path)
    artifacts = [manifest, _database_artifact(db_path, ".tsv")] + [
        Path(f"{base}.{part}.{extension}")
        for extension in ("bt2l", "bt2") for part in INDEX_PARTS
    ]
    if any(path.exists() for path in artifacts):
        parser.error(f"database artifacts already exist for {db_path}; choose a new -d path")
    if db_path.exists():
        parser.error(
            f"{db_path} exists without a manifest or index; if it is a failed earlier build, "
            "remove it, otherwise choose a new -d path"
        )


def _final_build_config(config: dict, db_path: Path) -> dict:
    final = config.copy()
    final["barcodemapper_db"] = str(db_path)
    final["barcodemapper_index_base"] = str(_database_base(db_path))
    final["barcodemapper_db_headers"] = str(_database_artifact(db_path, ".tsv"))
    return final


def _build_manual(args, parser: argparse.ArgumentParser, db_path: Path) -> None:
    _manual_target_guard(parser, db_path)
    missing = [tool for tool in ("snakemake", "bowtie2-build") if shutil.which(tool) is None]
    if missing:
        parser.error("required executables are missing from PATH: " + ", ".join(missing))

    bold = Path(args.bold).expanduser().resolve()
    unite = Path(args.unite).expanduser().resolve()
    stage_dir = db_path.parent / f".{db_path.name}.build"
    stage_path = stage_dir / "stage.json"
    tmp_dir = Path(args.tmpdir).expanduser().resolve() if args.tmpdir else stage_dir / "work"
    stage_key = {
        "bold_sha256": _sha256(bold),
        "unite_sha256": _sha256(unite),
        "bold_retain": args.bold_retain,
        "unite_retain": args.unite_retain,
        "version": __version__,
        "tmp_dir": str(tmp_dir),
    }
    reused = False
    if stage_dir.exists():
        try:
            previous_key = json.loads(stage_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous_key = None
        if previous_key != stage_key:
            shutil.rmtree(stage_dir)
        else:
            reused = True
    stage_dir.mkdir(parents=True, exist_ok=True)
    stage_path.write_text(json.dumps(stage_key, indent=2) + "\n", encoding="utf-8")

    staged_db = stage_dir / db_path.name
    tmp_dir.mkdir(parents=True, exist_ok=True)
    config = _build_config(staged_db, tmp_dir, args, bold=bold, unite=unite)
    config_path = tmp_dir / "data" / "config.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(config, sort_keys=True), encoding="utf-8")
    try:
        run_snakemake(tmp_dir, config_path, "slurm" if args.slurm else "local")
    except subprocess.CalledProcessError:
        if reused:
            unlock_command = ["barcodemapper", "--unlock", "-t", str(tmp_dir)]
            if args.slurm:
                unlock_command.append("--slurm")
            print(f"If a stale Snakemake lock blocks the resumed build, run: {shlex.join(unlock_command)}",
                  file=sys.stderr)
        raise

    incomplete = [path for path in _build_outputs(staged_db) if not path.is_file() or path.stat().st_size == 0]
    if incomplete:
        raise RuntimeError("database build is missing staged artifacts: " + ", ".join(map(str, incomplete)))
    _write_manifest(_database_artifact(staged_db, ".manifest.json"), staged_db,
                    _final_build_config(config, db_path), {}, {}, database_path=db_path,
                    source_info={
                        "bold": {"path": str(bold), "sha256": stage_key["bold_sha256"]},
                        "unite": {"path": str(unite), "sha256": stage_key["unite_sha256"]},
                    }, automatic_sources=False)
    _promote_database(staged_db, db_path)
    shutil.rmtree(stage_dir)
    print(f"Built database: {db_path}")
    print(f"Run manifest: {_database_artifact(db_path, '.manifest.json')}")


def _build_latest(args, parser: argparse.ArgumentParser, db_path: Path) -> None:
    try:
        releases = discover_latest()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(f"could not discover the latest BOLD/UNITE releases: {exc}")
    manifest = _auto_manifest(db_path)
    if db_path.exists() and manifest is None:
        parser.error(f"{db_path} has no automatic-build manifest; choose a new -d path")
    tmp_dir = Path(args.tmpdir).expanduser().resolve() if args.tmpdir else Path(
        "barcodemapper_" + datetime.now().strftime("%Y%m%d%H%M%S%f")
    ).resolve()
    same = _same_auto_database(db_path, manifest, releases, args)
    if same:
        index_files = [Path(f"{_database_base(db_path)}.{part}.bt2l") for part in INDEX_PARTS]
        headers = _database_artifact(db_path, ".tsv")
        index_current = _current_files(index_files, db_path)
        headers_current = _current_files([headers], db_path)
        if index_current and headers_current:
            print(f"Database is current: {db_path}")
            return
        missing = [tool for tool in ("snakemake", "bowtie2-build") if shutil.which(tool) is None]
        if missing:
            parser.error("required executables are missing from PATH: " + ", ".join(missing))
        tmp_dir.mkdir(parents=True, exist_ok=True)
        config = _build_config(db_path, tmp_dir, args, generate_index=not index_current,
                               generate_headers=not headers_current)
        config_path = tmp_dir / "data" / "config.yaml"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        config_path.write_text(yaml.safe_dump(config, sort_keys=True), encoding="utf-8")
        run_snakemake(tmp_dir, config_path, "slurm" if args.slurm else "local")
        print(f"Repaired database artifacts: {db_path}")
        return

    missing = [tool for tool in ("snakemake", "bowtie2-build") if shutil.which(tool) is None]
    if missing:
        parser.error("required executables are missing from PATH: " + ", ".join(missing))
    db_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{_database_base(db_path).name}.build-", dir=db_path.parent) as stage_name:
        stage_dir = Path(stage_name)
        staged_db = stage_dir / db_path.name
        try:
            paths, provenance = download_sources(releases, stage_dir / "sources")
        except (OSError, ValueError, json.JSONDecodeError, tarfile.TarError) as exc:
            parser.error(f"could not download the latest BOLD/UNITE FASTAs: {exc}")
        config = _build_config(staged_db, tmp_dir, args, bold=paths["bold"], unite=paths["unite"])
        config_path = tmp_dir / "data" / "config.yaml"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        config_path.write_text(yaml.safe_dump(config, sort_keys=True), encoding="utf-8")
        run_snakemake(tmp_dir, config_path, "slurm" if args.slurm else "local")
        manifest_config = _final_build_config(config, db_path)
        manifest_config.pop("bold_db")
        manifest_config.pop("unite_db")
        _write_manifest(_database_artifact(staged_db, ".manifest.json"), staged_db, manifest_config, {}, {},
                        source_info=provenance, database_path=db_path, automatic_sources=True)
        _promote_database(staged_db, db_path)
    print(f"Built latest database: {db_path}")
    print(f"Run manifest: {_database_artifact(db_path, '.manifest.json')}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Map shotgun reads to COI and ITS reference sequences.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-i", "--input", help="Directory of paired .fq.gz or .fastq.gz files")
    parser.add_argument("-1", "--read1", help="Comma-separated R1 FASTQ files")
    parser.add_argument("-2", "--read2", help="Comma-separated R2 FASTQ files")
    parser.add_argument("-o", "--output", help="Output taxonomy TSV file")
    parser.add_argument("-d", "--database", help="Combined BarcodeMapper FASTA path (.fa.gz for builds)")
    parser.add_argument("-b", "--bold", help="BOLD FASTA for database creation")
    parser.add_argument("-u", "--unite", help="UNITE FASTA for database creation")
    parser.add_argument("-x", "--bold_retain", default="k__Animalia", help="Comma-separated BOLD taxa to retain")
    parser.add_argument("-y", "--unite_retain", default="k__Viridiplantae,c__Agaricomycetes,c__Pezizomycetes", help="Comma-separated UNITE taxa to retain")
    parser.add_argument("-t", "--tmpdir", help="Working directory; reuse this path to resume")
    parser.add_argument("-s", "--slurm", action="store_true", help="Submit workflow jobs through Slurm")
    parser.add_argument("--unlock", action="store_true", help="Unlock an existing working directory")
    parser.add_argument("--build", action="store_true", help="Build and index a database only")
    parser.add_argument("-m", "--max_mismatches", "--max-mismatches", type=float, default=3.0,
                        help="Maximum NM edit distance as a percent of aligned bases (default: 3)")
    parser.add_argument("-c", "--min_coverage", "--min-coverage", type=int, default=80,
                        help="Minimum aligned bases (default: 80)")
    hits = parser.add_mutually_exclusive_group()
    hits.add_argument("-k", "--max-hits", type=int, default=DEFAULT_MAX_HITS,
                      help=f"Maximum alignments Bowtie2 reports per read (default: {DEFAULT_MAX_HITS})")
    hits.add_argument("--all-hits", "--full", dest="all_hits", action="store_true",
                      help="Report all Bowtie2 alignments per read (-a) instead of capping them; "
                           "slower and needs more memory")
    parser.add_argument("--no-bam", action="store_true",
                        help="Do not keep the Bowtie2 alignments as BAM files in <tmpdir>/bowtie (saves disk space)")
    parser.add_argument("--adapter-r1", default=DEFAULT_ADAPTER, help="R1 adapter sequence, or 'auto'")
    parser.add_argument("--adapter-r2", default=DEFAULT_ADAPTER, help="R2 adapter sequence, or 'auto'")
    parser.add_argument("--min-read-length", type=int, default=100, help="Minimum read length after fastp (default: 100)")
    parser.add_argument("--complexity-threshold", type=int, default=30, help="fastp complexity threshold (default: 30)")
    parser.add_argument("--no-low-complexity-filter", action="store_true", help="Disable fastp low complexity filtering")
    presets = parser.add_mutually_exclusive_group()
    for name, (percent, bases) in PRESETS.items():
        presets.add_argument(f"--{name}", action="store_true", help=f"Set -m {percent:g} and -c {bases}")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _parser()
    args = parser.parse_args(argv)
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        parser.exit(2)

    if args.build and args.unlock:
        parser.error("--build and --unlock cannot be combined")
    if args.build:
        if args.input or args.read1 or args.read2 or args.output:
            parser.error("--build accepts database inputs only; omit reads and -o")
    elif not args.unlock:
        if not args.output:
            parser.error("-o/--output is required for an analysis run")
        if bool(args.input) == bool(args.read1 or args.read2):
            parser.error("provide either -i or both -1 and -2")
        if bool(args.read1) != bool(args.read2):
            parser.error("both -1 and -2 are required for file-list input")
    if args.unlock and not args.tmpdir:
        parser.error("--unlock requires the existing --tmpdir")
    if args.unlock:
        if shutil.which("snakemake") is None:
            parser.error("snakemake is missing from PATH")
        tmp_dir = Path(args.tmpdir).expanduser().resolve()
        config_path = tmp_dir / "data" / "config.yaml"
        if not config_path.is_file():
            parser.error(f"no existing run configuration found at {config_path}")
        unlock_snakemake(tmp_dir, config_path, "slurm" if args.slurm else "local")
        return
    if not args.database:
        parser.error("-d/--database is required")
    if bool(args.bold) != bool(args.unite):
        parser.error("provide both BOLD and UNITE FASTA files, or neither")
    if not math.isfinite(args.max_mismatches) or args.max_mismatches < 0 or args.max_mismatches > 100:
        parser.error("-m/--max-mismatches must be between 0 and 100")
    if args.min_coverage < 1 or args.min_read_length < 1:
        parser.error("coverage and read length must be positive")
    if args.all_hits:
        args.max_hits = None
    if args.max_hits is not None and args.max_hits < 1:
        parser.error("-k/--max-hits must be positive")
    if not 0 <= args.complexity_threshold <= 100:
        parser.error("--complexity-threshold must be between 0 and 100")
    for label, adapter in (("--adapter-r1", args.adapter_r1), ("--adapter-r2", args.adapter_r2)):
        if adapter != "auto" and not re.fullmatch(r"[ACGTNacgtn]+", adapter):
            parser.error(f"{label} must be a DNA sequence or 'auto'")
    if (args.adapter_r1 == "auto") != (args.adapter_r2 == "auto"):
        parser.error("use 'auto' for both adapters or provide both sequences")

    db_path = Path(args.database).expanduser().resolve()
    if not db_path.name.lower().endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        parser.error("database path must end in .fa, .fasta, .fa.gz, or .fasta.gz")
    if args.build and db_path.suffix.lower() != ".gz":
        if db_path.exists():
            if args.bold:
                _manual_target_guard(parser, db_path)
            parser.error(f"uncompressed database target already exists: {db_path}; choose a new -d path")
        db_path = Path(f"{db_path}.gz")
    elif not args.build and not db_path.is_file() and db_path.suffix.lower() != ".gz":
        compressed_path = Path(f"{db_path}.gz")
        if compressed_path.is_file():
            db_path = compressed_path
    if args.build and db_path.suffix.lower() == ".gz":
        uncompressed_path = db_path.with_suffix("")
        if uncompressed_path.exists() and not db_path.exists():
            if args.bold:
                _manual_target_guard(parser, uncompressed_path)
            parser.error(f"uncompressed database already exists at {uncompressed_path}; choose a new -d basename")
    if args.build and not args.bold:
        _build_latest(args, parser, db_path)
        return
    if args.bold:
        for label, value in (("BOLD", args.bold), ("UNITE", args.unite)):
            if not Path(value).expanduser().is_file():
                parser.error(f"{label} FASTA does not exist: {value}")
        if args.build:
            db_path.parent.mkdir(parents=True, exist_ok=True)
            _build_manual(args, parser, db_path)
            return
        if db_path.exists():
            parser.error(f"database target already exists: {db_path}; choose a new -d path")
        db_path.parent.mkdir(parents=True, exist_ok=True)
    elif not db_path.is_file():
        parser.error(f"database FASTA does not exist: {db_path}")
    if not args.bold:
        try:
            validate_database(db_path)
        except (ValueError, OSError, EOFError) as exc:
            parser.error(f"invalid database FASTA: {exc}")

    tmp_dir = Path(args.tmpdir).expanduser().resolve() if args.tmpdir else Path(
        "barcodemapper_" + datetime.now().strftime("%Y%m%d%H%M%S%f")
    ).resolve()
    index_base, index_extension, generate_index, headers_path, generate_headers = _database_artifacts(
        db_path, tmp_dir, args.build
    )
    required_tools = ["snakemake"]
    if generate_index:
        required_tools.append("bowtie2-build")
    if not args.build:
        required_tools += ["bowtie2", "samtools", "fastp"]
    missing = [tool for tool in required_tools if shutil.which(tool) is None]
    if missing:
        parser.error("required executables are missing from PATH: " + ", ".join(missing))

    for preset, (percent, bases) in PRESETS.items():
        if getattr(args, preset):
            args.max_mismatches, args.min_coverage = percent, bases
            break

    tmp_dir.mkdir(parents=True, exist_ok=True)
    output_path = Path(args.output).expanduser().resolve() if args.output else None
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path = output_path.with_suffix(".html") if output_path else None
    config_path = tmp_dir / "data" / "config.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config = {
        "run": datetime.now(timezone.utc).isoformat(),
        "version": __version__,
        "barcodemapper_db": str(db_path),
        "barcodemapper_index_base": str(index_base),
        "barcodemapper_index_extension": index_extension,
        "generate_index": generate_index,
        "barcodemapper_db_headers": str(headers_path),
        "generate_headers": generate_headers,
        "tmp_dir": str(tmp_dir),
        "output_file": str(output_path) if output_path else "no",
        "build_only": "yes" if args.build else "no",
        "report_file": str(report_path) if report_path else "none",
        "max_edit_percent": args.max_mismatches,
        "min_coverage": args.min_coverage,
        "max_hits": args.max_hits,
        "keep_bam": not args.no_bam,
        "adapter_r1": args.adapter_r1,
        "adapter_r2": args.adapter_r2,
        "min_read_length": args.min_read_length,
        "complexity_threshold": args.complexity_threshold,
        "low_complexity_filter": not args.no_low_complexity_filter,
    }
    if args.bold:
        config.update({
            "bold_db": str(Path(args.bold).expanduser().resolve()),
            "unite_db": str(Path(args.unite).expanduser().resolve()),
            "bold_retain": args.bold_retain,
            "unite_retain": args.unite_retain,
        })

    read1: dict[str, list[str]] = {}
    read2: dict[str, list[str]] = {}
    if not args.build:
        try:
            if args.input:
                read1, read2 = inputdir_to_samples(args.input, tmp_dir)
            else:
                read1, read2 = inputlist_to_samples(args.read1, args.read2, tmp_dir)
        except ValueError as exc:
            parser.error(str(exc))
    config_path.write_text(yaml.safe_dump(config, sort_keys=True), encoding="utf-8")
    run_snakemake(tmp_dir, config_path, "slurm" if args.slurm else "local")
    manifest_path = (_database_artifact(db_path, ".manifest.json") if args.build
                     else output_path.with_suffix(".manifest.json"))
    _write_manifest(manifest_path, db_path, config, read1, read2)
    print(f"Run manifest: {manifest_path}")


if __name__ == "__main__":
    main()
