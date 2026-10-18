# Changelog

Versions follow [Semantic Versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`). The package version is declared in `pyproject.toml`, mirrored in `CITATION.cff`, and checked against this file before a release tag is created.

## [Unreleased]

### Added

- Add `-k` / `--max-hits N` to cap Bowtie2 alignments per read.
- Add `--all-hits` (alias `--full`) to report all Bowtie2 alignments (`-a`), as earlier versions did.

- Keep the Bowtie2 alignments as BAM files in `<tmpdir>/bowtie/`, written with `--no-unal` (unaligned reads are left out), so thresholds can be revisited without remapping. Use `--no-bam` to skip them.

### Changed

- Bowtie2 now reports at most 200 alignments per read by default (`-k 200`). Benchmarks showed 0-1.1% of read assignments changing, unchanged detected taxa, and a much faster, lighter mapping step. Use `--all-hits` for the previous behaviour; results can differ slightly from 1.1.x otherwise.
- The mapping step requests more memory (6x instead of 3x the input size) and time (3x) under `--all-hits`, since uncapped runs ran out of memory or time at the capped request.

## [1.1.3] - 2026-10-16

This release also ships 1.1.1 and 1.1.2, which were never tagged or published.

### Fixed

- Read mapping now requests resources sized from a measured run instead of the previous guess: 3 times the input size in memory (peak use was about 2.3 times) and 12 minutes per GB of input (about 3 hours with the default database), doubled on each retry. The old 5 minutes per GB (about 74 minutes) made every 100,000-pair run against the default database time out on its first attempt, and the 5 times memory request was about twice what mapping used once `samtools collate` was removed.

### Changed

- Read mapping no longer pipes Bowtie2 through `samtools collate`. With a large reference (for example 10.5 million sequences, a 1.2 GB `@SQ` header) collate copied the full header into each of its 64 temporary shards, so mapping needed about 100 GB of memory and 15 GB of temporary files even for 10,000 read pairs. The alignments now stream straight to BAM with `samtools view -b`, and the filtered extraction table is sorted by read ID in the extract rule so each fragment's alignments stay adjacent for taxonomy assignment. Results are unchanged.

## [1.1.2] - 2026-10-11

### Fixed

- Manual `--build` with `-b` and `-u` now builds in a hidden staging directory beside the target and promotes the FASTA, Bowtie2 index, headers, and manifest only when all are complete, so a failed build no longer leaves a partial database.
- A failed manual build keeps its stage and resumes on rerun with the same inputs and work directory; the stage is discarded if the BOLD/UNITE checksums, retain options, version, or work directory change. A stale Snakemake lock can be cleared with `barcodemapper --unlock -t <work directory>`.
- Builds into a target that already holds database artifacts stop with a clear message to choose a new `-d` path; a lone FASTA asks for review before removal or reuse of that path.

## [1.1.1] - 2026-10-10

### Added

- Install instructions for the release wheel without a source checkout.
- A Bioconda recipe with the workflow executables as runtime dependencies.
- Reproducible Mjolnir benchmark plans for mock communities, reference novelty, method comparisons, and scaling.

### Fixed

- Slurm workflows no longer cap rule threads at the controller's CPU count; the Slurm profile now sets `cores: 8`, so `index_db` and read mapping run with their intended threads.
- `index_db` memory and runtime are sized from the uncompressed database rather than the gzipped file, so large reference builds no longer exhaust their retries with out-of-memory and timeout failures.

## [1.1.0] - 2026-10-10

The first tagged release. Version 1.0.0 was an untagged development snapshot; its changes are included here.

### Added

- Snakemake workflow for database preparation, read filtering, Bowtie2 mapping, taxonomy assignment, aggregation, and HTML reporting.
- Local and Slurm profiles, mapping presets, and paired FASTQ input support.
- Package metadata, MIT license, software citation, release checks, and build CI.
- End-to-end synthetic workflow tests and a run manifest with tool versions, parameters, database and source checksums, and sample filtering counts.
- Documentation for installation, source provenance, limitations, validation, and release steps.
- Allow `barcodemapper --build` without source files to fetch the current BOLD and UNITE releases, record their provenance, and skip or repair an already current database.

### Changed

- Map with Bowtie2's exhaustive alignment mode to preserve ambiguous taxonomic matches.
- Count edit distance against actual aligned bases and expose adapter, read-length, and complexity controls.
- Keep conflicting taxonomies for identical BOLD sequences so lowest-common-ancestor assignment can represent ambiguity.
- Make the HTML report self-contained and include zero-hit and kingdom-only assignments in its plots.
- Validate paired reads, numeric inputs, and database paths before starting Snakemake.
- Compress newly built reference FASTAs with gzip and read compressed databases directly for validation, indexing, and taxonomy overviews.

### Fixed

- Use fastp's low-complexity flag with its threshold argument correctly.
- Preserve empty samples and unresolved taxonomy during aggregation and reporting.
- Pass paths safely to workflow commands and propagate Bowtie2 pipe failures.
- Generate missing analysis index and taxonomy overview files in the work directory instead of beside a shared reference FASTA.
- Handle sample names that match internal report column names.
