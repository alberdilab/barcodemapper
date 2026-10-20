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

Each level filters alignments independently after one Bowtie2 mapping pass. `-c` / `--min-coverage` sets the minimum aligned bases per alignment; `-m` / `--max-mismatches` sets the maximum Bowtie2 `NM` edit distance as a percentage of those aligned bases. `NM` includes substitutions and indels. By default, all four presets run in this order:

| Preset | Minimum aligned bases | Maximum edits |
| --- | ---: | ---: |
| `--ultrastrict` | 120 | 0% |
| `--strict` | 100 | 2% |
| `--standard` | 80 | 3% |
| `--permissive` | 60 | 5% |

Use `--stringency` with a comma-separated list to select presets, or `--stringency all` for all four. Preset flags add levels to that selection. With no `--stringency` or preset flags, `-m` and/or `-c` define one `custom` level instead of the four defaults; with presets, they add `custom` as an extra level. An omitted custom threshold defaults to 3% or 80 aligned bases, respectively. Multiple levels are reported in preset order, followed by `custom`.

```sh
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv --stringency strict,standard
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv --strict --permissive
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv -m 1 -c 110
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv --strict -m 1 -c 110
```

The first command writes `results/taxonomy.strict.tsv` and `results/taxonomy.standard.tsv`; the second writes `results/taxonomy.strict.tsv` and `results/taxonomy.permissive.tsv`. The third writes one custom table to `results/taxonomy.tsv`; the fourth writes `results/taxonomy.strict.tsv` and `results/taxonomy.custom.tsv`. A single preset level, such as `--strict` alone, also writes directly to `-o`. Ultrastrict requires at least 120 aligned bases, while the default minimum read length after fastp is 100; some retained reads therefore cannot pass ultrastrict. Changing which alignments are accepted can change the lowest common ancestor, so taxonomy counts at each rank are not simple subsets across levels.

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

Common analysis options are listed below.

| Option | Purpose |
| --- | --- |
| `-i`, `--input` | Directory of paired `.fq.gz` / `.fastq.gz` files |
| `-1`, `--read1` / `-2`, `--read2` | Comma-separated R1 / R2 files |
| `-o`, `--output` | Taxonomy TSV path for one level; base path for level-named TSVs when several levels run |
| `-d`, `--database` | Combined BarcodeMapper FASTA (default: newest `bmdb_*` in `--db-dir` or the current directory) |
| `--db-dir`, `--db-name` | Where to find databases, and which one (folder name or date) |
| `-b`, `--bold` / `-u`, `--unite` | Source FASTAs for `--build` |
| `-x`, `--bold_retain` / `-y`, `--unite_retain` | Taxa retained when building |
| `-t`, `--tmpdir` | Working directory; reuse to resume |
| `-s`, `--slurm` | Submit jobs through Slurm |
| `--build`, `-n`, `--keep` | Build and index a database only; name it; prune older ones |
| `--list`, `--verify`, `--update`, `--check` | Inspect, check or refresh databases (see [](database.md)) |
| `--unlock` | Unlock a working directory |
| `--stringency`, preset flags, `-m`, `-c` | Select presets or add a custom stringency level |
| `-k`, `--all-hits`, `--no-bam` | Alignment cap, uncapped mapping, BAM retention |
| `--version` | Print the version |
