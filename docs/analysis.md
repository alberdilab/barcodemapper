# Running an analysis

Put one read pair per sample in a directory, named `sample_1.fq.gz` and `sample_2.fq.gz` (or `.fastq.gz`):

```sh
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv
```

## Inputs

Instead of `-i`, provide comma-separated R1 and R2 file lists with `-1` and `-2`. Both mates must exist and use matching sample names. Merge lanes before analysis, or name each pair as a separate sample. The names `Taxonomy`, `Phylum`, `Kingdom`, and `color` are reserved for output columns.

## Working directory and resuming

`-t work_directory` gives the run a stable working directory that Snakemake can resume. The default is a new timestamped directory. Use `--unlock -t work_directory` only if a previous run left a stale Snakemake lock.

## Read preprocessing

Reads are filtered with fastp. The bundled adapter is `CTGTCTCTTATACACATCT` for both mates. Set `--adapter-r1` and `--adapter-r2` for your library, or set both to `auto` to use fastp's adapter handling.

| Option | Default | Meaning |
| --- | --- | --- |
| `--adapter-r1`, `--adapter-r2` | `CTGTCTCTTATACACATCT` | Adapter sequence, or `auto` |
| `--min-read-length` | 100 | Minimum read length after fastp |
| `--complexity-threshold` | 30 | fastp complexity threshold |
| `--no-low-complexity-filter` | off | Disable low-complexity filtering |

Inspect fastp's per-sample HTML and JSON files in the work directory for filtering results.

## Mapping stringency

`-c` / `--min-coverage` is the minimum aligned bases per alignment; `-m` / `--max-mismatches` is the maximum Bowtie2 `NM` edit distance as a percentage of those aligned bases. `NM` includes substitutions and indels. The default is 80 bases and 3%. Presets replace both values:

| Preset | Minimum aligned bases | Maximum edits |
| --- | ---: | ---: |
| `--ultrastrict` | 120 | 0% |
| `--strict` | 100 | 2% |
| `--standard` (default) | 80 | 3% |
| `--permissive` | 60 | 5% |

```sh
barcodemapper -i reads --db-dir references -o results.tsv --strict
```

## Alignment cap (`-k`) and `--all-hits`

By default Bowtie2 reports at most 200 alignments per read (`-k 200`). In benchmarks on six samples this changed 0-1.1% of read assignments against reporting all hits, left the detected taxa unchanged, and cut the mapping step from up to five hours to a few minutes and peak memory by about half.

- `-k` / `--max-hits N` changes the cap. Bowtie2 `-k` stops after finding N valid alignments, which are not necessarily the best N, so for reads with more than N equally good hits the shared taxonomy can be assigned from a subset.
- `--all-hits` (alias `--full`) reports every alignment (`-a`). This can be much slower on repetitive sequences and large databases, and the pipeline requests more memory (6x the input size instead of 3x) and time (3x) for the mapping step automatically.
- `--all-hits` cannot be combined with `-k`.

## Keeping alignments as BAM

Bowtie2 alignments are kept as one BAM per sample in `<tmpdir>/bowtie/`, written with `--no-unal` so reads that did not align are left out and the files stay small. Because `-m` and `-c` are applied after Bowtie2, each BAM holds every alignment Bowtie2 reported, including those a stricter threshold excluded, so the data can be re-extracted at other thresholds without remapping.

Alignments are still limited by Bowtie2's own minimum score and by the `-k` cap, and the BAMs do not record how many reads failed to align. Use `--no-bam` to skip them and save disk space.

## Slurm

Use `--slurm` on a cluster with the Snakemake Slurm executor plugin installed. Adjust `barcodemapper/profile/slurm/config.yaml` for site resources before running large datasets; the provided profile submits up to 50 jobs and caps threads per job at 8.

## Option reference

Run `barcodemapper --help` for the authoritative list.

| Option | Purpose |
| --- | --- |
| `-i`, `--input` | Directory of paired `.fq.gz` / `.fastq.gz` files |
| `-1`, `--read1` / `-2`, `--read2` | Comma-separated R1 / R2 files |
| `-o`, `--output` | Output taxonomy TSV |
| `-d`, `--database` | Combined BarcodeMapper FASTA (default: newest `bmdb_*` in `--db-dir` or the current directory) |
| `--db-dir`, `--db-name` | Where to find databases, and which one (folder name or date) |
| `-b`, `--bold` / `-u`, `--unite` | Source FASTAs for `--build` |
| `-x`, `--bold_retain` / `-y`, `--unite_retain` | Taxa retained when building |
| `-t`, `--tmpdir` | Working directory; reuse to resume |
| `-s`, `--slurm` | Submit jobs through Slurm |
| `--build`, `-n`, `--keep` | Build and index a database only; name it; prune older ones |
| `--list`, `--verify`, `--update`, `--check` | Inspect, check or refresh databases (see [](database.md)) |
| `--unlock` | Unlock a working directory |
| `-m`, `-c`, presets | Mapping stringency |
| `-k`, `--all-hits`, `--no-bam` | Alignment cap, uncapped mapping, BAM retention |
| `--version` | Print the version |
