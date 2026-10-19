# Mjolnir execution record

## 2026-10-10 startup

The benchmark workspace is
`/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks`.
The v1.1.0 wheel was transferred and verified against `wheel.sha256` on Mjolnir
before job submission. The release tag resolves to commit
`7efbdd4941a74ac3b7b5e24ca420be6ba075aed1`.

| Job | Stage | Dependency |
| --- | --- | --- |
| 45898208 | Create Project 01 Conda environment and explicit package lock | none |
| 45898209 | Freeze BOLD and UNITE source FASTAs and SHA-256 hashes | 45898208 |
| 45898213 | Download and verify fungal shotgun run `SRR14351851` | none |
| 45898216 | Build and index the default reference | 45898209 |
| 45898217 | Build and index the all-fungi reference | 45898209 |
| 45898219 | Create a synchronized 100,000-pair fungal pilot | 45898213 |
| 45898220 | Run BarcodeMapper `--standard` on the fungal pilot | 45898217 and 45898219 |
| 45898237 | Hash and freeze the default reference artifacts | 45898216 |
| 45898238 | Hash and freeze the all-fungi reference artifacts | 45898217 |
| 45898244 | Clone Project 02 environment from the explicit lock | none |
| 45898245 | Clone Project 03 environment from the explicit lock | 45898244 |
| 45898246 | Clone Project 04 environment from the explicit lock | 45898244 |
| 45898268 | Download and verify animal mock `ERR6001161` | none |
| 45898269 | Download and verify probable extraction blank `ERR6001173` | none |
| 45898270 | Download and verify fungal pooled biomass `SRR14351849` | none |
| 45898272 | Create 100,000-pair animal mock pilot | 45898268 |
| 45898273 | Create 100,000-pair extraction-control pilot | 45898269 |
| 45898274 | Run the default-reference animal mock pilot | 45898216 and 45898272 |
| 45898275 | Run the default-reference extraction-control pilot | 45898216 and 45898273 |
| 45898284 | Verify and record the fungal pilot | 45898220 and 45898238 |
| 45898285 | Verify and record the animal mock pilot | 45898274 and 45898237 |
| 45898286 | Verify and record the extraction-control pilot | 45898275 and 45898237 |
| 45898392 | Run Project 04's 100,000-pair real-shotgun scaling pilot | 45898274 and 45898237 |
| 45898413 | Inventory default-reference taxonomy | 45898237 |
| 45898414 | Inventory all-fungi reference taxonomy | 45898238 |
| 45898560 | Sample Project 04 pilot scratch occupancy every 30 seconds | after 45898392 starts |
| 45898599 | Initial Project 02 source-index attempt; found empty BOLD sequences | none |
| 45898603 | Index nonempty source accessions for Project 02 holdout design | 45898209 |
| 45898609 | Make checksum-traced source copies excluding empty FASTA records | 45898209 |
| 45898620–45898633 | Replacement reference builds, audits, pilots, and scaling pilot | see [recovery job map](metadata/references/job-map.tsv) |
| 45898641 | Annotate expected taxa against verified reference inventories | 45898624 and 45898625 |
| 45898643 | Project 03 best-accepted-hit ablation on the fungal pilot | 45898627 |
| 45898644 | Project 03 best-accepted-hit ablation on the animal pilot | 45898629 |

The FASTQ download uses Mjolnir's `filetransfer` partition and QoS. The source
snapshot performs decompression and uses `cpuqueue`. The build and pilot
controllers use `cpuqueue` and submit their
Snakemake workers through the packaged Slurm profile. Slurm submission commands
and output logs are under `tasks/` in the workspace. Check current status with
`squeue -u jpl786` and final states with `sacct -j JOBID`.

The 2026-10-10 ENA run reports, their retrieval URLs and hashes, sample XML,
and inclusion decisions are retained in [metadata/01-mock-communities](metadata/01-mock-communities/).
The selection rule is paired Illumina WGS. It includes all 12 animal mock
replicates and `MC_C1`, predeclared as the probable extraction blank based on
the alias, smaller run, and the paper's one-blank design. The fungal runs are
`SRR14351851` (pooled DNA) and `SRR14351849` (pooled biomass). Both papers'
Table 1 member names, sample roles, source URLs, and XML hashes were frozen
before any classification result was available.
All 15 selected ENA runs have download jobs. The final 11 animal mock transfers
are serialized to limit concurrent traffic; their run-to-job mapping is
[`remaining-jobs.tsv`](metadata/01-mock-communities/remaining-jobs.tsv).

Environment setup job `45898208` completed successfully (exit code 0). Its
Linux package lock, package JSON, YAML and hashes are retained under
`metadata/01-mock-communities/software/`.
Environment clones `45898244`, `45898245`, and `45898246` completed
successfully. All four explicit locks have the same SHA-256
`6cd8bf42a7a2597547f8bfe343ea6a2b730d80455b44f7fe7e98ff6106bf03d8`.
The YAML, locks, package JSON and checksums for each project are in `metadata/`.

Source snapshot job `45898209` completed successfully. The frozen source record
and checksum manifest are under `metadata/references/`. They identify BOLD
Public `09-Oct-2026` and UNITE All Eukaryotes `2025-02-19`; the uncompressed
source FASTAs are 16,981,910,548 and 202,353,058 bytes respectively. The
source files on Mjolnir are read-only. Both reference build controllers have
started and submitted their first worker jobs. They were cancelled after a
source defect was confirmed, as described below; no reference snapshot was
completed by those jobs.
The build controllers initially logged `sacct` connection warnings. Their
compute-node Slurm configuration pointed accounting to `localhost`; the login
node's configuration successfully queries the cluster accounting service. A
project-local copy of that configuration and [sacct-wrapper.sh](sacct-wrapper.sh)
now provide the working setting for `sacct` without changing the pinned wheel
or its Snakemake profile. The wrapper is linked into all four Conda prefixes.
Its configuration and checksums are retained in `metadata/slurm/`.
The controller warning logs stopped growing after the wrapper was installed;
confirm controller completion when the workers finish.

The BOLD snapshot includes **2,336,565 zero-length records** among 24,491,027
records. A small v1.1.0 reproduction confirmed that its database preparation
keeps these records and its validator then rejects the reference. To avoid a
known failure after a long build, I cancelled the two original build
controllers and workers (45898216, 45898217, 45898292, 45898293), their
dependent reference and pilot jobs, and the dependent scaling pilot. The
cancelled workers' partial `.fa.gz` outputs were removed after confirming
that no complete manifest or index existed; their logs remain under `tasks/`.
The [source filter](drop_empty_source_records.py) retains the original files
and copies every nonempty record byte-for-byte into immutable derivatives.
Job `45898609` completed successfully. Its
[manifest](metadata/references/empty-record-filter.json) records original and
derived SHA-256 hashes, record counts, and the exclusion-table hash. UNITE had
zero empty records. The complete excluded-header table remains on Mjolnir
because it is 261 MB. The replacement job chain is in
[job-map.tsv](metadata/references/job-map.tsv), with builds starting after the
filter and downstream pilots waiting on verified reference snapshots.

Project 02's source-index retry `45898603` completed successfully. The
[index summary](metadata/02-reference-novelty/source-index.json) records
20,408,974 BOLD COI rows and 168,030 fungal UNITE ITS rows before
length/taxonomy eligibility filtering, plus the SQLite index checksum. The
7.8 GB index remains on Mjolnir. Project 01's first animal mock download and
synchronized 100,000-pair pilot input also completed, with their manifests
retained in `metadata/01-mock-communities/`. Once both reference audits finish,
job `45898641` will write per-run expected binomials and exact-name reference
coverage before classification outputs are assessed. The fungal paper has
43 expected rows after excluding its possible contaminant, but only 38
distinct binomial names because several strains share species names.
Project 03's [best-hit ablation](best_hit_ablation.py) is queued for both
verified pilot datasets. It uses the same fastp-filtered reads and accepted
alignment TSV as BarcodeMapper, selecting lowest `NM`, then longest aligned
coverage, with a lexical accession tie-break. The output includes one row per
retained pair, including no-hit pairs, and records source and output hashes.

At 2026-10-10 08:11 UTC, both replacement reference controllers (`45898620`,
`45898621`) and the next animal transfer (`45898512`) were running. The
reference finalizers, audits, pilot analyses, expected-taxon annotation,
best-hit ablations, and scaling pilot were waiting on their dependencies.
No reference build or classification result has yet been verified. The next
gate is successful reference finalization, followed by pilot output and
resource review before full-scale runs.

## 2026-10-10 BarcodeMapper 1.1.1 and reference build resubmission

Both reference builds (`45898620`, `45898621`) failed in `index_db`. The
Snakemake workers ran out of memory three times and then hit the 4 h time limit
on the final retry. Two defects in v1.1.0 caused this: its Slurm profile set no
`cores`, so Snakemake capped every rule at the controller's 2 hardware threads
(`bowtie2-build --threads 2` instead of 8), and `index_db` sized memory and
runtime from the gzipped FASTA (about 1.1 GB for roughly 7.9 Gbp uncompressed).
The prepared `default.fa.gz` and `all-fungi.fa.gz` were complete and are reused;
database preparation is unchanged between 1.1.0 and 1.1.1.

The fix is released as BarcodeMapper **1.1.1**. Its wheel
`barcodemapper-1.1.1-py3-none-any.whl` has SHA-256
`67f31651097153ba7310ae26e08c3486876635e303d5c19b7b7386d56c967a67`
(see `wheel.sha256`). It was built from the working tree on top of commit
`893040ff24eb85a706ce2b4bf97d12e55619a11d`, before the 1.1.1 commit existed.
It was installed with `pip install --no-deps --force-reinstall` into all four
project prefixes. The Conda locks are unchanged. The 1.1.0 wheels, checksums and
earlier plan files are retained under `superseded-1.1.0/` in each location.
No run produced by v1.1.0 is used in any result; the pilot inputs and cleaned
sources are unaffected.

The earlier chain was cancelled (`45898622`-`45898633`, `45898641`,
`45898643`, `45898644`) and its map kept as
`tasks/clean-source-recovery/job-map.superseded-1.1.0.tsv`. The replacement
chain is:

| Job | Stage | Dependency |
| --- | --- | --- |
| 45898822 | Build default reference | 45898609 |
| 45898823 | Build all-fungi reference | 45898609 |
| 45898824 | Freeze default reference | 45898822 |
| 45898825 | Freeze all-fungi reference | 45898823 |
| 45898826 | Audit default reference | 45898824 |
| 45898827 | Audit all-fungi reference | 45898825 |
| 45898828 | Fungal pilot | 45898825 |
| 45898829 | Verify fungal pilot | 45898828 |
| 45898830 | Animal pilot | 45898824 and 45898272 |
| 45898831 | Verify animal pilot | 45898830 |
| 45898832 | Blank pilot | 45898824 |
| 45898833 | Verify blank pilot | 45898832 |
| 45898834 | Scaling pilot, 100k pairs | 45898831 |
| 45898835 | Scratch sampler | 45898834 |
| 45898836 | Annotate expected taxa | 45898826 and 45898827 |
| 45898837 | Best-hit ablation, fungal | 45898829 |
| 45898838 | Best-hit ablation, animal | 45898831 |

The `index_db` memory and runtime constants are estimates from the failed run.
Check each build's peak memory and elapsed time with `sacct` once it finishes.

## 2026-10-11 BarcodeMapper 1.1.2 resubmission

Builds `45898822` and `45898823` of the 1.1.1 chain exited after 11 s with
"database target already exists": the prepared `default.fa.gz` and
`all-fungi.fa.gz` from the failed 1.1.1 build were still at the final path.
BarcodeMapper **1.1.2** stages manual builds in `<dir>/.<db-filename>.build/`,
resumes on rerun and promotes the database only when complete.

The wheel `barcodemapper-1.1.2-py3-none-any.whl` has SHA-256
`f67b4abea4384639f401db84c7d0f29e9bfa55984bc8b976b85c9a4649197c44`
(see `wheel.sha256`). It was installed with
`pip install --no-deps --force-reinstall` into all four project prefixes
(`barcodemapper --version` prints 1.1.2 in each). The 1.1.1 wheels, checksums
and plan checksum are retained under `superseded-1.1.1/` in each location.

The dead 1.1.1 chain (`45898824`-`45898838`) was cancelled; its map is
`tasks/clean-source-recovery/job-map.superseded-1.1.1.tsv`. Nothing was
deleted. Moved to `references/superseded-1.1.1/`: `default.fa.gz`,
`all-fungi.fa.gz`, their `.tsv` manifests, the partial `*.bt2l.tmp` index files,
and `work-default-clean`, `work-all-fungi-clean`. `references/cleaned-sources/`
is untouched. Everything is submitted by
[submit_clean_source_recovery_1.1.2.sh](submit_clean_source_recovery_1.1.2.sh);
the builds have no dependency because source cleaning (`45898609`) completed.

| Job | Stage | Dependency |
| --- | --- | --- |
| 45899205 | Build default reference | none |
| 45899206 | Build all-fungi reference | none |
| 45899207 | Freeze default reference | 45899205 |
| 45899208 | Freeze all-fungi reference | 45899206 |
| 45899209 | Audit default reference | 45899207 |
| 45899210 | Audit all-fungi reference | 45899208 |
| 45899211 | Fungal pilot | 45899208 |
| 45899212 | Verify fungal pilot | 45899211 |
| 45899213 | Animal pilot | 45899207 and 45898272 |
| 45899214 | Verify animal pilot | 45899213 |
| 45899215 | Blank pilot | 45899207 |
| 45899216 | Verify blank pilot | 45899215 |
| 45899217 | Annotate expected taxa | 45899209 and 45899210 |
| 45899218 | Best-hit ablation, fungal | 45899212 |
| 45899219 | Best-hit ablation, animal | 45899214 |
| 45899220 | Scaling pilot, 100k pairs | 45899214 |
| 45899221 | Scratch sampler | 45899220 (`after`) |

## Project 05 (not yet submitted)

Scripts: [05_taxid_paths.sbatch](05_taxid_paths.sbatch) builds both taxonomy
tables; [05_blast_nt.sbatch](05_blast_nt.sbatch) searches nt;
[05_kraken2_nt.sbatch](05_kraken2_nt.sbatch) classifies against Kraken2 nt;
[submit_05.sh](submit_05.sh) prints the dependency chain unless
`CONFIRM_SUBMIT=yes` is set. Jobs require `RUN` for BLAST and Kraken2.
The Kraken2 path build takes kingdom taxids from
`/maps/datasets/globe_databases/ncbi/taxonomy/20241001/nodes.dmp` because its
`ktaxonomy.tsv` uses numeric D/K codes for these nodes. Override that source
with `KRAKEN_KINGDOM_TAXDUMP` if needed; the audit records its checksum and
how many kingdom taxids occur in the Kraken2 taxonomy.
Optional variables are `READ_IDS` (BLAST subset), `THREADS`, `OUT_DIR`,
`PRELOAD` (Kraken2, default 1), and `MEM` (Kraken2 submission, default 1150G).
All scripts also accept path overrides (`BENCH_ROOT`, `BENCH_ENV`, `PLAN_DIR`,
`PROJECT_DIR`, database paths, taxonomy paths, and FASTQ paths). The submit
helper accepts `READ_IDS_<RUN>` for each pilot.

Copy the entire contents of `benchmarks/mjolnir`, including
`metadata/05-default-databases/name-map.tsv`, into `$BENCH_ROOT/plans/` before
submission. The two mock pilots require `READ_IDS_<RUN>` subset files. Review
the `--report` outputs from `05_taxid_paths.sbatch` and compare the name map
against `ktaxonomy.tsv` before reading Kraken2 results. A failed Kraken2 job
leaves later Kraken2 jobs pending with `DependencyNeverSatisfied` by design.
The Kraken2 preload is best effort: classifier-only timing includes any first
read of `hash.k2d` that the cache did not retain, so report the first run with
that caveat.
