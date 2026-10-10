![BarcodeMapper logo](https://raw.githubusercontent.com/alberdilab/barcodemapper/main/images/barcodemapper.png)

# BarcodeMapper

BarcodeMapper maps paired shotgun reads to taxonomically annotated COI and ITS reference sequences. It filters reads with fastp, aligns them with Bowtie2, and assigns each read ID to the lowest taxonomy shared by its accepted alignments. It produces a cumulative taxonomy count table, a self-contained HTML report, and a JSON run manifest.

Read counts describe marker-matching fragments. They are **not** estimates of organism abundance: marker copy number, reference coverage, read length, and mapping ambiguity affect the counts. See [validation and limitations](docs/validation.md).

## Install

Python 3.12 or newer is required. The provided Conda environment pins the workflow tools and Python dependencies. From a checked-out release tag or source tree:

```sh
conda env create --file barcodemapper_env.yml
conda activate barcodemapper
python -m pip install --no-deps .
barcodemapper --version
```

The pinned environment is intended for Linux x86_64 and macOS x86_64. The pinned fastp 0.23.1 build is unavailable for native macOS arm64; Apple Silicon users can use an x86_64 Conda environment under Rosetta or install compatible workflow executables separately. For a Python-only install, `python -m pip install .` installs the Python dependencies; `snakemake`, `fastp`, `bowtie2`, `bowtie2-build`, and `samtools` must also be on `PATH`. Slurm runs additionally need the Snakemake Slurm executor plugin.

## Reference database

BarcodeMapper requires a FASTA whose headers have the form `accession|marker|k__Kingdom;p__Phylum;...`, where `marker` is `COI-5P` or `ITS`. The FASTA and Bowtie2 index are separate from the Python package.

An existing March 2025 database snapshot was previously distributed as [barcodemapper_db_202503.tar.gz](https://sid.erda.dk/share_redirect/FGathhqKb5/barcodemapper_db_202503.tar.gz). Its source inputs were BOLD Public (14 March 2025) and UNITE All Eukaryotes (19 February 2025). Confirm its availability and the underlying providers' terms before using or redistributing it. The software's MIT license does not grant rights to third-party reference sequences.

To build your own database, obtain the FASTA inputs from [BOLD](https://www.boldsystems.org/) and the [UNITE repository](https://unite.ut.ee/repository.php), record their release identifiers, and run:

```sh
barcodemapper --build \
  -b BOLD_Public.fa -u UNITE_All_Eukaryotes.fasta \
  -d references/barcodemapper.fa \
  -x k__Animalia \
  -y k__Viridiplantae,c__Agaricomycetes,c__Pezizomycetes
```

`-x` and `-y` accept comma-separated taxonomy terms present in the source data. The default filters select Animalia from BOLD and Viridiplantae, Agaricomycetes, and Pezizomycetes from UNITE. The build creates the combined FASTA, six Bowtie2 index files, a taxonomy overview TSV, and a manifest with source and output SHA-256 hashes. It checks reference headers and requires at least one retained record for each marker. Use a new `-d` path for each build; existing FASTA files are not overwritten. Analysis runs also check the FASTA header format before mapping.

When reporting results, cite both the software release and the reference providers. Follow [BOLD's citation guidance](https://boldsystems.org/about/citation/) for the platform and any data package DOI, and [UNITE's citation guidance](https://unite.ut.ee/cite.php) for its current database paper. Retain the exact source release identifiers and checksums with your analysis.

## Run an analysis

Put one read pair per sample in a directory, named `sample_1.fq.gz` and `sample_2.fq.gz` (or `.fastq.gz`):

```sh
barcodemapper -i reads -d references/barcodemapper.fa -o results/taxonomy.tsv
```

You can instead provide comma-separated R1 and R2 file lists with `-1` and `-2`. Both mates must exist and use matching sample names. Merge lanes before analysis, or name each pair as a separate sample. The names `Taxonomy`, `Phylum`, `Kingdom`, and `color` are reserved for output columns. `-t work_directory` gives the run a stable working directory that Snakemake can resume. The default is a new timestamped directory. Use `--unlock -t work_directory` only if a previous run left a stale Snakemake lock.

The bundled adapters are `CTGTCTCTTATACACATCT` for both mates. Set `--adapter-r1` and `--adapter-r2` for your library, or set both to `auto` to use fastp's adapter handling. Other preprocessing controls include `--min-read-length`, `--complexity-threshold`, and `--no-low-complexity-filter`. Inspect fastp's per-sample HTML and JSON files in the work directory for filtering results.

### Mapping stringency

`-c` / `--min-coverage` is the minimum aligned bases per alignment; `-m` / `--max-mismatches` is the maximum Bowtie2 `NM` edit distance as a percentage of those aligned bases. `NM` includes substitutions and indels. The default is 80 bases and 3%. Presets replace both values:

| Preset | Minimum aligned bases | Maximum edits |
| --- | ---: | ---: |
| `--ultrastrict` | 120 | 0% |
| `--strict` | 100 | 2% |
| `--standard` (default) | 80 | 3% |
| `--permissive` | 60 | 5% |

For example, `barcodemapper -i reads -d references/barcodemapper.fa -o results.tsv --strict`. Bowtie2 reports all alignments it finds, so repetitive sequences and large databases can increase run time and disk use substantially.

### Slurm

Use `--slurm` on a cluster with the Snakemake Slurm executor plugin installed. Adjust `barcodemapper/profile/slurm/config.yaml` for site resources before running large datasets; the provided profile submits up to 50 jobs and caps threads per job at 8.

## Outputs

For `-o results/taxonomy.tsv`, BarcodeMapper writes:

| File | Contents |
| --- | --- |
| `results/taxonomy.tsv` | Taxonomy paths and one count column per sample, including zero-hit samples. |
| `results/taxonomy.html` | Offline, self-contained bar and sunburst plots. |
| `results/taxonomy.manifest.json` | Version, parameters, tool versions, database checksum, inputs, and fastp/assignment counts. |

Counts are propagated to every resolved ancestor. A fragment assigned to species therefore also appears in its genus, family, and higher-rank rows. **Do not sum all taxonomy rows** to obtain a total. A fragment with no shared named taxon is `Unclassified`; a sample with no accepted alignments has a zero column. The report's relative bars use assigned fragments at the phylum level plus unresolved and unclassified fragments as the denominator, and are labelled as a fraction of marker-matching fragments.

## Tests and releases

Run `python -m pytest -q` after installing `pytest` and the workflow tools. The two end-to-end tests skip when their external executables are unavailable. [Test details](tests.md) and [validation scope](docs/validation.md) explain what they cover.

The package version is declared once in `pyproject.toml` and exposed through `barcodemapper --version`. Versions follow `MAJOR.MINOR.PATCH`; release tags use `vMAJOR.MINOR.PATCH`. Changes are recorded in [CHANGELOG.md](CHANGELOG.md), citation metadata in [CITATION.cff](CITATION.cff), and release steps in [docs/releasing.md](docs/releasing.md). Version 1.0.0 is currently an untagged source snapshot; no package publication or release tag is implied.

BarcodeMapper source code is available under the [MIT license](LICENSE). Reference databases retain their providers' terms.
