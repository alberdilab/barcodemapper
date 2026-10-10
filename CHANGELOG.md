# Changelog

Versions follow [Semantic Versioning](https://semver.org/) (`MAJOR.MINOR.PATCH`). The package version is declared in `pyproject.toml`, mirrored in `CITATION.cff`, and checked against this file before a release tag is created. No release tag has been created yet.

## [Unreleased]

## [1.0.0] - Pending release

The first release is being prepared. This section describes the current source tree; no tag or package publication is implied.

### Added

- Snakemake workflow for database preparation, read filtering, Bowtie2 mapping, taxonomy assignment, aggregation, and HTML reporting.
- Local and Slurm profiles, mapping presets, and paired FASTQ input support.
- Modern package metadata, MIT license, software citation, release checks, and build CI.
- End-to-end synthetic workflow tests and a run manifest with tool versions, parameters, database and source checksums, and sample filtering counts.
- Documentation for installation, source provenance, limitations, validation, and release steps.

### Changed

- Map with Bowtie2's exhaustive alignment mode to preserve ambiguous taxonomic matches.
- Count edit distance against actual aligned bases and expose adapter, read-length, and complexity controls.
- Keep conflicting taxonomies for identical BOLD sequences so lowest-common-ancestor assignment can represent ambiguity.
- Make the HTML report self-contained and include zero-hit and kingdom-only assignments in its plots.
- Validate paired reads, numeric inputs, and database paths before starting Snakemake.

### Fixed

- Use fastp's low-complexity flag with its threshold argument correctly.
- Preserve empty samples and unresolved taxonomy during aggregation and reporting.
- Pass paths safely to workflow commands and propagate Bowtie2 pipe failures.
