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
| NCBI taxonomy | `ncbi/taxonomy/20251202` | Newer than both databases. Check taxon IDs in the Kraken2 `ktaxonomy.tsv` before relying on it. |

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
uses its own `ktaxonomy.tsv` and the BLAST arm uses `taxonomy/20251202`, and the
name-mapping table is built per arm. Do not translate Kraken2 IDs with the newer
dump.

**Kraken2 memory:** a 953 GB hash needs about 1 TB of RAM to load fully. The
2 TB nodes (for example `mjolnircomp08fl`, `mjolnircomp09fl`, 96 cores) can hold
it, so request one with at least 1.1 TB and measure peak memory. Do not use
`--memory-mapping` for timing runs unless the file is cached, because the first
read from the network file system would dominate the time. Stage or pre-read
the file and state which was done.

## Taxonomy harmonisation

NCBI names and BarcodeMapper's BOLD/UNITE paths differ, so the common schema
in [convert_calls.py](convert_calls.py) cannot be reused as is.

- Freeze an NCBI-to-common name mapping **before** looking at any accuracy
  result. Resolve synonyms with the NCBI `names.dmp` and record every mapping.
- Evaluate only at ranks where a name can be mapped without ambiguity
  (expected: phylum, class, order, family, genus; species where unambiguous).
  Report how many mock taxa are lost to the mapping at each rank.
- Add `nt-blastn` and `nt-kraken2` subcommands (and `nt-centrifuge` if built)
  to `convert_calls.py`, with tests on small fixtures that cover synonyms,
  unmapped names, and unclassified reads.
- Calls that map to non-target organisms (bacteria, host, human) are counted
  as off-target calls, not as errors at the marker level.

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

## BLAST subset (decided)

Run BLAST nt on all three pilots. Use the same predeclared stratified subset and
seeds as Project 03 for the two mock pilots, so the nt and BOLD/UNITE results
compare read for read. For the blank, search every fastp-retained pair if the
count is small (cap it at the Project 03 subset size otherwise), because
false-positive rate is the point of that control and sampling would hide rare
calls. Record the subset files and their hashes.

Centrifuge is not part of this project.

## Deliverables and completion check

Save the nt provenance record, the name-mapping table, build commands, converted
calls, accuracy and blank-control tables, and timing logs. The project is
complete when every arm used the Project 03 FASTQ hashes, the database provenance
is recorded, and the mapping table was frozen before the accuracy results.
