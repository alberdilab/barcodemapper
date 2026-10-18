# Reference database

BarcodeMapper requires a FASTA whose headers have the form

```text
accession|marker|k__Kingdom;p__Phylum;...
```

where `marker` is `COI-5P` or `ITS`. The FASTA and Bowtie2 index are separate from the Python package. Database builds write a gzip-compressed FASTA (`.fa.gz` or `.fasta.gz`); analysis also accepts existing uncompressed FASTAs.

An existing March 2025 snapshot was previously distributed as [barcodemapper_db_202503.tar.gz](https://sid.erda.dk/share_redirect/FGathhqKb5/barcodemapper_db_202503.tar.gz). Its source inputs were BOLD Public (14 March 2025) and UNITE All Eukaryotes (19 February 2025). Confirm its availability and the underlying providers' terms before using or redistributing it. The software's MIT license does not grant rights to third-party reference sequences.

## Database folders

Databases live in dated folders. By default a build creates `bmdb_<YYYYMMDD>/` in the current directory, holding everything the database needs:

```text
bmdb_20261018/
├── bmdb_20261018.fa.gz          # combined FASTA
├── bmdb_20261018.*.bt2l         # Bowtie2 index (six files)
├── bmdb_20261018.tsv            # taxonomy overview
└── bmdb_20261018.manifest.json  # sources, filters, checksums
```

| Option | Effect |
| --- | --- |
| `--db-dir DIR` | Build: the folder in which `bmdb_<YYYYMMDD>/` is created. Mapping, `--list`, `--verify`: where to look. `DIR` may hold several `bmdb_*` folders or be one of them (`/projects/refs` or `/projects/refs/bmdb_20261018/`). Default: the current directory. |
| `-n, --name NAME` | Label the database: `bmdb_<YYYYMMDD>_<NAME>`. Letters, digits, `.`, `_` and `-` only. |
| `-d, --database FILE` | Use an explicit FASTA path instead (for builds, `.fa.gz`). Cannot be combined with `--db-dir`. |

## Automatic build

Download the current BOLD Public snapshot and the current UNITE All Eukaryotes general FASTA release (including singletons), then build and index the database:

```sh
barcodemapper --build                       # ./bmdb_<YYYYMMDD>/
barcodemapper --build --db-dir /projects/refs -n coi-its
barcodemapper --build -d references/barcodemapper.fa.gz   # explicit path
```

This mode checks the providers' current release identifiers every time it runs:

- If the FASTA checksum, source releases, filters, and index/overview files are all current, it exits without downloading or rebuilding.
- If only an index or overview file is missing or stale, it regenerates that artifact.
- If a source release or filter has changed, it downloads the new sources, builds in a temporary directory, and replaces the previous automatic build only after the new build succeeds.

It will not replace an existing database that has no automatic-build manifest; use a new `-d` path in that case. Source release identifiers, DOIs, download URLs, and FASTA SHA-256 hashes are recorded in the manifest. The BOLD snapshot is several gigabytes, so allow enough disk space for the extracted source FASTAs, the compressed combined FASTA, and the Bowtie2 index. Provider access rules and availability can change.

## Finding the database during analysis

Without `-d`, an analysis run uses the newest `bmdb_*` database in `--db-dir` (default: the launch directory). If none is found it stops and asks for `-d` or `--db-dir`.

```sh
barcodemapper -i reads -o results/taxonomy.tsv                      # newest in ./
barcodemapper -i reads -o results/taxonomy.tsv --db-dir /projects/refs
barcodemapper -i reads -o results/taxonomy.tsv --db-dir /projects/refs --db-name 20260101
```

`--db-name` pins an older database by folder name or date (`bmdb_20260101`, `20260101`), which is useful for reproducing an earlier analysis. The run manifest records the database folder name, build date and source releases alongside the checksum.

## Managing databases

```sh
barcodemapper --list --db-dir /projects/refs      # date, sources, filters, index, size
barcodemapper --verify --db-dir /projects/refs    # newest; add --db-name to pick one
barcodemapper --update --check --db-dir /projects/refs
barcodemapper --update --db-dir /projects/refs --keep 3
```

- `--list` summarises every database from its manifest.
- `--verify` checks the FASTA headers, the Bowtie2 index, the overview file and the manifest checksum, and exits with status 1 if any check fails.
- `--update` queries BOLD and UNITE and builds a new dated database only when a source release or filter differs from the newest local one (`-n NAME` compares against databases with that name). With `--check` it only reports, exiting with status 10 when a newer release exists, so it can run from a scheduler.
- `--keep N` (with `--build` or `--update`, without `-d`) deletes all but the newest `N` databases of the same name in `--db-dir` after a successful build. Only folders that contain a manifest are removed.

## Build from your own FASTA files

To use files you obtained yourself from [BOLD](https://www.boldsystems.org/) and the [UNITE repository](https://unite.ut.ee/repository.php):

```sh
barcodemapper --build \
  -b BOLD_Public.fa -u UNITE_All_Eukaryotes.fasta \
  -d references/barcodemapper.fa.gz \
  -x k__Animalia \
  -y k__Viridiplantae,c__Agaricomycetes,c__Pezizomycetes
```

`-x` and `-y` accept comma-separated taxonomy terms present in the source data. The defaults select Animalia from BOLD and Viridiplantae, Agaricomycetes, and Pezizomycetes from UNITE.

The build creates:

- the combined gzip FASTA (compression level 9);
- six Bowtie2 index files (Bowtie2 indexes the compressed FASTA directly);
- a taxonomy overview TSV;
- a manifest with source and output SHA-256 hashes.

If `-d` ends in `.fa` or `.fasta` for a build, `.gz` is appended. The build checks reference headers and requires at least one retained record for each marker. For builds from user-supplied FASTAs, use a new `-d` path each time: a target that already holds a database is refused, and a lone FASTA without an index or manifest must be reviewed before removal or reuse of that path.

### Staged, resumable builds

Manual builds work in a hidden directory `.<database filename>.build/` next to the target (for example `references/.barcodemapper.fa.gz.build/`), and `--tmpdir` defaults to its `work/` subdirectory. The FASTA, Bowtie2 index, headers overview, and manifest are moved into place only when all are complete, and the stage is then removed.

If a build fails, the stage is kept; rerun the same command to resume. The stage is keyed on the BOLD and UNITE SHA-256 hashes, the retain options, the BarcodeMapper version, and the work directory, and is discarded and rebuilt if any of them change. If a stale Snakemake lock blocks a resumed build:

```sh
barcodemapper --unlock -t references/.barcodemapper.fa.gz.build/work
```

Add `--slurm` if the build used it, or give the same `-t` path used for the build.

## Index reuse during analysis

Analysis checks the FASTA header format before mapping. It reuses a complete, current Bowtie2 index (`.bt2l` or `.bt2`) and taxonomy overview beside the FASTA when available. If either is missing or older than the FASTA, it generates that artifact under the run's working directory instead, so analysis never creates files beside a shared reference FASTA. Reuse `-t work_directory` to reuse an index generated there on a resumed run; separate work directories each need their own index unless the database has a prebuilt one.

## Citing the references

When reporting results, cite both the software release and the reference providers. Follow [BOLD's citation guidance](https://boldsystems.org/about/citation/) for the platform and any data package DOI, and [UNITE's citation guidance](https://unite.ut.ee/cite.php) for its current database paper. Retain the exact source release identifiers and checksums with your analysis.
