![BarcodeMapper logo](https://raw.githubusercontent.com/alberdilab/barcodemapper/main/images/barcodemapper.png)

# BarcodeMapper

BarcodeMapper maps paired shotgun reads to taxonomically annotated COI and ITS reference sequences. It filters reads with fastp, aligns them with Bowtie2, and assigns each read ID to the lowest taxonomy shared by its accepted alignments. It produces a cumulative taxonomy count table, a self-contained HTML report, and a JSON run manifest.

Read counts describe marker-matching fragments. They are **not** estimates of organism abundance: marker copy number, reference coverage, read length, and mapping ambiguity affect the counts. See [validation and limitations](docs/validation.md).

## Install

Python 3.12 or newer is required. To install the released wheel without cloning the repository, create an environment with the workflow tools and install the wheel:

```sh
conda create -n barcodemapper --override-channels -c conda-forge -c bioconda \
  python=3.12 snakemake=8.27.1 fastp bowtie2 samtools \
  pandas pyyaml biopython pysam plotly pip
conda activate barcodemapper
python -m pip install --no-deps \
  https://github.com/alberdilab/barcodemapper/releases/download/v1.1.0/barcodemapper-1.1.0-py3-none-any.whl
barcodemapper --version
```

The provided Conda environment pins all workflow tools and Python dependencies. From a checked-out release tag or source tree, it can also be installed with:

```sh
conda env create --file barcodemapper_env.yml
conda activate barcodemapper
python -m pip install --no-deps .
barcodemapper --version
```

The pinned environment is intended for Linux x86_64 and macOS x86_64. The pinned fastp 0.23.1 build is unavailable for native macOS arm64; Apple Silicon users can use an x86_64 Conda environment under Rosetta or install compatible workflow executables separately. For a Python-only install, `python -m pip install .` installs the Python dependencies; `snakemake`, `fastp`, `bowtie2`, and `samtools` must also be on `PATH`. Database builds and analyses without a current Bowtie2 index also need `bowtie2-build` (included with Bowtie2). Slurm runs additionally need the Snakemake Slurm executor plugin. The release wheel above is installable with pip now; the short `pip install barcodemapper` and `conda install -c bioconda barcodemapper` commands require publication to PyPI and Bioconda, respectively.

## Reference database

BarcodeMapper requires a FASTA whose headers have the form `accession|marker|k__Kingdom;p__Phylum;...`, where `marker` is `COI-5P` or `ITS`. The FASTA and Bowtie2 index are separate from the Python package. Database builds write a gzip-compressed FASTA (`.fa.gz` or `.fasta.gz`); analysis also accepts existing uncompressed FASTAs.

An existing March 2025 database snapshot was previously distributed as [barcodemapper_db_202503.tar.gz](https://sid.erda.dk/share_redirect/FGathhqKb5/barcodemapper_db_202503.tar.gz). Its source inputs were BOLD Public (14 March 2025) and UNITE All Eukaryotes (19 February 2025). Confirm its availability and the underlying providers' terms before using or redistributing it. The software's MIT license does not grant rights to third-party reference sequences.

To download the current BOLD Public snapshot and the current UNITE All Eukaryotes general FASTA release (including singletons), then build and index the database:

```sh
barcodemapper --build -d references/barcodemapper.fa.gz
```

This mode checks the providers' current release identifiers every time it runs. If the FASTA checksum, source releases, filters, and index/overview files are already current, it exits without downloading or rebuilding. If only an index or overview file is missing or stale, it regenerates that artifact. If a source release or filter has changed, it downloads the new sources, builds in a temporary directory, and replaces the previous automatic build only after the new build succeeds. It will not replace an existing database without an automatic-build manifest; use a new `-d` path in that case. Source release identifiers, DOIs, download URLs, and FASTA SHA-256 hashes are recorded in the manifest. The BOLD snapshot is several gigabytes, so allow enough disk space for the extracted source FASTAs, compressed combined FASTA, and Bowtie2 index during the build. Provider access rules and availability can change.

To use FASTA files you obtained yourself from [BOLD](https://www.boldsystems.org/) and the [UNITE repository](https://unite.ut.ee/repository.php), run:

```sh
barcodemapper --build \
  -b BOLD_Public.fa -u UNITE_All_Eukaryotes.fasta \
  -d references/barcodemapper.fa.gz \
  -x k__Animalia \
  -y k__Viridiplantae,c__Agaricomycetes,c__Pezizomycetes
```

`-x` and `-y` accept comma-separated taxonomy terms present in the source data. The default filters select Animalia from BOLD and Viridiplantae, Agaricomycetes, and Pezizomycetes from UNITE. The build creates the combined gzip FASTA at compression level 9, six Bowtie2 index files, a taxonomy overview TSV, and a manifest with source and output SHA-256 hashes. Bowtie2 indexes the compressed FASTA directly. If `-d` ends in `.fa` or `.fasta` for a build, `.gz` is appended to the output path. It checks reference headers and requires at least one retained record for each marker. For builds from user-supplied FASTAs, use a new `-d` path each time; existing FASTA files are not overwritten. Analysis runs also check the FASTA header format before mapping.

An analysis reuses a complete, current Bowtie2 index (`.bt2l` or `.bt2`) and taxonomy overview beside the FASTA when available. If either is missing or older than the FASTA, it generates that artifact under the run's working directory instead. Analysis therefore does not create files beside a shared reference FASTA. Reuse `-t work_directory` to reuse an index generated there on a resumed run; separate work directories require their own index unless the database has a prebuilt one.

When reporting results, cite both the software release and the reference providers. Follow [BOLD's citation guidance](https://boldsystems.org/about/citation/) for the platform and any data package DOI, and [UNITE's citation guidance](https://unite.ut.ee/cite.php) for its current database paper. Retain the exact source release identifiers and checksums with your analysis.

## Run an analysis

Put one read pair per sample in a directory, named `sample_1.fq.gz` and `sample_2.fq.gz` (or `.fastq.gz`):

```sh
barcodemapper -i reads -d references/barcodemapper.fa.gz -o results/taxonomy.tsv
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

For example, `barcodemapper -i reads -d references/barcodemapper.fa.gz -o results.tsv --strict`. Bowtie2 reports all alignments it finds, so repetitive sequences and large databases can increase run time and disk use substantially.

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

The package version is declared once in `pyproject.toml` and exposed through `barcodemapper --version`. Versions follow `MAJOR.MINOR.PATCH`; release tags use `vMAJOR.MINOR.PATCH`. Changes are recorded in [CHANGELOG.md](CHANGELOG.md), citation metadata in [CITATION.cff](CITATION.cff), and release steps in [docs/releasing.md](docs/releasing.md). Version 1.1.0 is the first tagged release; 1.0.0 was an untagged development snapshot.

BarcodeMapper source code is available under the [MIT license](LICENSE). Reference databases retain their providers' terms.
