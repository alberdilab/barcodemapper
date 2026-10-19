# Project 05: practitioner defaults on NCBI nt

## Question

How does BarcodeMapper, with its BOLD and UNITE reference, compare with what a
user would get from general-purpose shotgun classifiers run against NCBI nt?
[Project 03](03-method-comparison.md) holds the reference constant to isolate
the algorithm. This project changes the reference as well, so it measures the
whole practitioner workflow. Never pool its results with Project 03 and never
read a difference here as an algorithm effect.

## Arms

| Arm | Database | Notes |
| --- | --- | --- |
| BarcodeMapper v1.2.3 | Frozen BOLD and UNITE reference | Use the 1.2.3 release and confirm `barcodemapper --version` prints 1.2.3. Reuse the Project 03 1.2.3 pilot rerun if it exists; otherwise run the three pilots once here. Do not pool with 1.1.2 results. |
| BLASTn 2.16.0 | Shared NCBI `nt` BLAST database, `ncbi/nt/20250722` | Same frozen settings as [blast_config.yaml](blast_config.yaml), on the Project 03 subset. |
| Kraken2 2.1.3 | Shared prebuilt `kraken2/kraken2_nt/20240530` | Primary k-mer arm. Not rebuilt. See the note on memory below. |

Kraken2 with its Standard database (RefSeq bacteria, archaea, viruses, human)
is not an arm. It holds almost no COI or ITS, so it only shows that the
default fails. Run it once on the pilots as a one-line sanity check and report
the unclassified fraction, nothing more.

## Inputs

- The same fastp-filtered paired FASTQs, with the same hashes, as Project 03:
  fungal SRR14351851, animal ERR6001161 and the blank ERR6001173.
- Project 01 taxon truth for the real mocks. Project 02 per-fragment truth is
  not available, so there is no synthetic-truth evaluation here.
- Run COI and ITS reads separately, as in Project 03.

## Database provenance

Use the existing shared databases under `/maps/datasets/globe_databases/` and
do not download a new nt. Checked on Mjolnir:

| Use | Path | Notes |
| --- | --- | --- |
| BLAST nt | `ncbi/nt/20250722` | Newest full nt. Also available: `ncbi/nt_core/20251119`. |
| Kraken2 nt | `kraken2/kraken2_nt/20240530` | Newest; `hash.k2d` is 953 GB. The 20231129 build is 1.5 TB with its tarball. |
| NCBI taxonomy | `ncbi/taxonomy/20241001` (decided 2026-10-19) | Used for the BLAST arm and as the kingdom source for the Kraken2 path map. It matches the Kraken2 database (May 2024) and BarcodeMapper's reference classes. `20251202` shows the same taxa under newer classes and orders (for example Pichiomycetes, Serinales) and is kept only as an optional sensitivity table. |

1. Treat these as read-only. Copy nothing; point to them.
2. Record each path, directory date, file list and the `nt.md5` / `*.md5`
   checksums, and the Kraken2 `opts.k2d` (k-mer and minimizer settings).
3. The three databases come from different dates (Kraken2 May 2024, BLAST July
   2025, taxonomy December 2025). Report the dates with every result and use
   each tool's own matching taxonomy, so the arms are not on one taxonomy.
4. The BOLD snapshot (`09-Oct-2026`) is newer than all of them. Some BOLD
   sequences are also in GenBank, so nt and BOLD overlap, but nt may lack
   records BarcodeMapper has. Note this when interpreting recall.
5. Read access was confirmed on 2026-10-18 as `jpl786`: `opts.k2d`,
   `ktaxonomy.tsv`, the checksum files, the BLAST `nt.000.nin`/`.nsq` listing and
   the `taxonomy/20251202` dumps all read normally despite the restrictive
   permission bits (an ACL grants access). The multi-hundred-GB files were not
   read.

**Taxonomy check (2026-10-18):** of 2,056,473 taxon IDs in the Kraken2
`ktaxonomy.tsv`, 9,764 (0.5%) are absent from `taxonomy/20251202/nodes.dmp`:
9,340 were merged into other IDs and 421 were deleted. Decision: the Kraken2 arm
uses its own `ktaxonomy.tsv` and the BLAST arm uses `taxonomy/20241001` (the 20251202 check above was done before this decision), and the
name-mapping table is built per arm. Do not translate Kraken2 IDs with the newer
dump.

**Kraken2 memory:** a 953 GB hash needs about 1 TB of RAM to load fully. The
2 TB nodes (for example `mjolnircomp08fl`, `mjolnircomp09fl`, 96 cores) can hold
it, so request one with at least 1.1 TB and measure peak memory. Do not use
`--memory-mapping` for timing runs unless the file is cached, because the first
read from the network file system would dominate the time. Stage or pre-read
the file and state which was done.

## Taxonomy harmonisation

NCBI names and BarcodeMapper's BOLD/UNITE paths differ, so each arm's calls are
mapped to the common `k__;p__;c__;o__;f__;g__;s__` form by
[build_ncbi_path_map.py](build_ncbi_path_map.py) (taxid to path, used by the
`kraken2` subcommand) and [build_nt_accession_paths.py](build_nt_accession_paths.py)
(BLAST accession to path, used by `blastn`). `convert_calls.py` is reused as is.
Calls at taxa with no domain ancestor (root, "cellular organisms", "synthetic
construct") become `unclassified`.

Survey of the reference against `taxonomy/20251202` (2026-10-18), distinct names
matching an NCBI scientific name at the same rank (default / all-fungi):
kingdom 2 of 3 and 1 of 2 (see below), phylum 34 of 48 / 38 of 50, class
109 of 151 / 143 of 203, order 666 of 793 / 748 of 924, family 5,847 of 6,571 /
5,959 of 6,794, genus 66,638 of 77,559 / 61,976 of 72,975. About 166,000
species labels are BOLD BIN placeholders that NCBI cannot match, which is
expected. Kingdoms: BOLD animals are `Animalia` (NCBI scientific name
`Metazoa`; `Animalia` is only a synonym), plants are `Viridiplantae` (matches
NCBI), fungi are `Fungi` (matches).

Rules, frozen before looking at any accuracy result:

- The default [name-map](metadata/05-default-databases/name-map.tsv) has one
  row, `Metazoa` to `Animalia`. Plants are **not** renamed.
- The k slot is the domain name for Bacteria, Archaea and Viruses; the NCBI
  `kingdom` node only applies under Eukaryota. Recent dumps put `kingdom` nodes
  (for example Bacillati, Pseudomonadati) under Bacteria.
- Many reference names differ from NCBI at class and order (for example
  Teleostei is an infraclass, Anthophyta and Eudicotyledonae do not exist).
  Same-rank matching would score these as wrong calls, so extend the name map
  with a reviewed entry for **each expected mock taxon** at those ranks, using
  the Project 01 `expected-members.tsv`. Build the list by joining the expected
  taxa to NCBI by name across ranks and synonyms, review it by hand, and freeze
  it with its SHA-256 in the run record.
- Report, for every arm, how many expected taxa could not be mapped at each
  rank. Evaluate a rank only where the mapping is unambiguous; unmatched called
  names count as wrong at that rank and are listed.
- Run the same comparison against the Kraken2 `ktaxonomy.tsv` before the
  Kraken2 arm, because that arm uses its own taxonomy.

The converters need tests; they are in `tests/test_build_ncbi_path_map.py`,
`tests/test_build_nt_accession_paths.py` and `tests/test_nt_integration.py`.

**Name-map review (2026-10-19), expected mock taxa.** All 66 expected species
exist in NCBI. Animal run (26 compared): only the kingdom differs systematically
(`Animalia` vs `Metazoa`, mapped); 3 species differ at family because they are
different taxa, and `Parachela_hn` vs `Parachela` at order is mapped (the
reference holds only the `_hn` form). Fungal run (33 compared): 18 differ at class
and 20 at order against `20251202`, because BarcodeMapper's reference and the
Kraken2 `ktaxonomy.tsv` use the older Saccharomycetes/Saccharomycetales
classification. This is version drift, not naming, and is the reason for choosing
`20241001`. Seven expected species have no exact BarcodeMapper species path
(spelling variants or renamed species) and are reported as such, not mapped.
The Kraken2 `ktaxonomy.tsv` codes kingdoms as `D1`/`D2`, so the kingdom taxids
come from the taxdump's `nodes.dmp` (`--kingdom-from-taxdump`).

## Reference and input facts (Mjolnir, 2026-10-18)

`B` is `/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks`.

- The fungal pilot SRR14351851 used `B/references/all-fungi.fa.gz`; the animal
  ERR6001161 and blank ERR6001173 used `B/references/default.fa.gz`. Record
  which reference each pilot uses in every run record.
- fastp-filtered pairs are in
  `B/projects/01-mock-communities/work/<RUN>-pilot-standard/fastp/<RUN>_{1,2}.fq.gz`.
  They come from the 1.1.2 environment (fastp 0.23.1). Use these files, hash
  them once, and give the same files to every arm and version.
- The 1.2.3 reference rebuild was still running at the time of the survey, so
  the 1.2.3 pilots cannot start until it is finished and audited.
- Environment with blastn 2.16.0 and kraken2 2.1.3:
  `B/projects/03-method-comparison/software/env-1.2.3`.
- The 2 TB nodes `mjolnircomp08fl` and `mjolnircomp09fl` (partition `cpuqueue`)
  were nearly fully allocated, so a Kraken2 nt job will queue.
- The shared databases are on NFS (Isilon). Nothing is cached on compute nodes,
  so the Kraken2 hash is read from the network on first use.

## Evaluation

1. Real mocks: taxa detected, unexpected calls, and fragment support under the
   Project 01 thresholds. Keep reference-absent mock members visible in the
   denominator table, and separately note members that are absent from nt.
2. Blank control: report calls per million pairs for each arm. This is the
   false-positive comparison nt is most likely to affect.
3. Compute: classifier-only wall time, CPU-hours, peak memory, index size and
   build cost, using [time_classifier.sh](time_classifier.sh) with the
   Project 03 rules. The nt build is much larger than the BOLD/UNITE build, so
   report it separately.
4. Use [evaluate_calls.py](evaluate_calls.py) with `--expected-taxa`, after
   confirming the harmonised paths load without changes. If it needs changes,
   add tests first.

## BLAST subset (decided, built 2026-10-19)

Project 02 truth does not exist and the pilots are almost entirely BarcodeMapper
no-hit (fungal 92,203 of 92,383 pairs; animal 98,134 of 98,155; blank 98,876 of
98,880), so the equal-strata sampler [make_blast_subset.py](make_blast_subset.py)
cannot be used. [make_pilot_blast_subset.py](make_pilot_blast_subset.py) builds
each subset from the 1.1.2 per-fragment taxonomy: **every** pair BarcodeMapper
classified, plus 5,000 no-hit pairs drawn per seed (seeds 1337, 2027, 4099).
BLAST nt runs once on the union of the seeds (14,394 / 14,271 / 14,257 pairs for
fungal / animal / blank), and each seed is evaluated as a subset of the union.
The blank is capped the same way, not searched in full.

Strata come from BarcodeMapper's own calls, not from truth. The no-hit draw is
what exposes calls BLAST makes that BarcodeMapper does not. Weight results by
stratum size when estimating overall rates, and say so in the write-up. The
lists and `manifest.json` (input and output hashes) are under
`projects/05-default-databases/inputs/blast-subsets/<RUN>/` on Mjolnir.

Kraken2 runs on all retained pairs, not on the subset. Centrifuge is not part of
this project.

## Deliverables and completion check

Save the nt provenance record, the name-mapping table, build commands, converted
calls, accuracy and blank-control tables, and timing logs. The project is
complete when every arm used the Project 03 FASTQ hashes, the database provenance
is recorded, and the mapping table was frozen before the accuracy results.
