# Project 03: controlled method comparison

## Question

What does BarcodeMapper's all-accepted-alignments and lowest-common-ancestor
assignment contribute compared with a single best hit and a k-mer classifier?
Separate classification accuracy from the cost of each search strategy.

## Inputs and fairness rules

Use the frozen FASTQs and taxon truth from [Project 01](01-mock-communities.md)
and the locked per-fragment truth from [Project 02](02-reference-novelty.md).
Use the **same reference records and taxonomy tree** for each method within
each dataset. Run COI and ITS separately in the read-level comparison and
compare mock-community taxa only at ranks represented in each database.

The first comparison uses the three verified pilots: fungal SRR14351851, animal
ERR6001161, and the blank ERR6001173, which is a negative control. Project 02
truth is not built yet, so the synthetic-truth evaluation is deferred until it
is. Best-hit outputs exist for the first two pilots only; the blank has none.
BarcodeMapper v1.2.3 is the primary version throughout, at its released
default (`-k 200` alignments per read). The earlier 1.1.2 pilot is history
only and does not appear in the comparison tables.

The default and all-fungi references are rebuilt with 1.2.3 from the retained
cleaned sources. Their FASTA hashes must match the 1.1.2 ones before use. fastp,
Bowtie2 and samtools are the versions pinned by 1.2.3 (see the 1.2.3 entry in
[CHANGELOG.md](../../CHANGELOG.md)), installed in a new Project 03 environment
with the wheel built from tag v1.2.3. The wheel SHA-256 is recorded in the run
record. Kraken2, BLAST, vsearch and minimap2 are installed in that environment.

All methods receive the same fastp-filtered paired FASTQs. Generate these
once with BarcodeMapper's pinned fastp settings and save their hashes; the
files under its `work/fastp/` directory are suitable once the primary run has
finished. Record both:

- **Classifier-only cost:** mapping/classification after fastp, excluding
  download, database build, fastp, and report rendering.
- **End-to-end cost:** raw FASTQ to final table, including preprocessing and
  report steps. State which timing is used in every plot.

## Methods

| Method | Preparation | Main output |
| --- | --- | --- |
| BarcodeMapper v1.2.3 (primary) | Wheel built from tag v1.2.3 (SHA-256 recorded in the run record); confirm that `barcodemapper --version` prints 1.2.3. Released default, `-k 200` per read, on the 1.2.3-rebuilt reference | Per-fragment taxonomy from `work/taxonomy/` |
| BarcodeMapper v1.2.3 `--all-hits` (secondary) | Same wheel and reference, pilots only, with `--all-hits`. Report separately and never pool with the primary row | Per-fragment taxonomy; shows the cost of the `-k 200` cap |
| Single best accepted Bowtie2 hit | From BarcodeMapper's `work/extract/` TSV, retain the alignment with the highest prespecified score; resolve exact ties with a fixed accession sort rule | Per-fragment taxonomy; isolates the assignment rule without changing mapping |
| Kraken2 2.1.3 | Build a custom database from the *same* FASTA records with deterministic synthetic taxon IDs for BarcodeMapper's taxonomy paths | Paired classification converted back to the common taxonomy paths |
| vsearch | `--usearch_global` against the same marker FASTA; keep the top hit per mate and take the LCA of tied top hits | Per-fragment taxonomy from a global-identity search |
| minimap2 (`-x sr`) | Index the same marker FASTA; keep all reported hits per pair and take their LCA | Per-fragment taxonomy from a second all-hits mapper |
| BLASTn 2.16.0, secondary | Index the same marker FASTA; search a predeclared stratified subset, including positives, off-targets, and ambiguous reads | Alignment cross-check, not the primary throughput comparator |

The best-hit ablation needs a deterministic score definition. BarcodeMapper's
extract TSV contains `NM` and covered bases, but no Bowtie2 alignment score;
either use Bowtie2's `AS` tag from its name-collated BAM, or define and label a
proxy such as lowest `NM` followed by longest coverage. Predeclare this choice
and tie handling before examining accuracy. Do not call this ablation a faster
pipeline: it reuses BarcodeMapper's all-alignments BAM. If a speed comparison
to Bowtie2 default reporting is desired, run that mapping mode separately.
For the queued pilot ablation, [best_hit_ablation.py](best_hit_ablation.py)
predeclares the proxy as lowest `NM`, then longest covered bases; exact score
ties use lexical accession, marker, and taxonomy order. It joins against every
fastp-retained R1 read ID so no-hit pairs remain in the output.

Build the Kraken2 inputs with [build_kraken_inputs.py](build_kraken_inputs.py)
`<barcodemapper-reference.fa.gz> <output-dir>`. It writes `taxonomy/nodes.dmp`,
`taxonomy/names.dmp`, `taxid-to-path.tsv`, `library.fa`, `seqid2taxid.map` and
`audit.json`, and fails unless per-marker record counts match the source FASTA.
It refuses to overwrite existing files. [submit_kraken_build.sbatch](submit_kraken_build.sbatch)
runs it, then `kraken2-build`, and records the commands in `build-commands.sh`.
It needs `BENCH_ROOT`, `BENCH_ENV`, `REFERENCE_FASTA`, `KRAKEN_DB`, `K`,
`MINIMIZER_LEN` and `MAX_DB_SIZE` exported before submission.

Kraken2 accepts [custom taxonomy assignments](https://kraken2.readthedocs.io/en/latest/custom_databases.html).
Create one taxon ID per distinct path prefix, with stable IDs from a sorted
list. Give every reference record its leaf ID without deduplicating conflicting
sequence labels. Archive the ID-to-path map and the complete Kraken2 database
build command. Audit that input sequence counts and per-marker totals exactly
match the BarcodeMapper FASTA. Use Kraken2 `--paired` and retain individual
read-pair classifications. A differing taxonomic ontology would confound the
comparison; this synthetic tree avoids using unrelated NCBI names or synonyms.

For BLASTn, define its task, e-value, identity, coverage, hit reporting, and
LCA/tie policy in a frozen configuration. Bowtie2 and BLASTn have different
search and scoring policies, so do not imply that equal numeric cutoffs mean
identical search sensitivity. Limit BLASTn to a predeclared subset large enough
to cover every challenge stratum and at least three seeds. The frozen settings
are in [blast_config.yaml](blast_config.yaml). Draw the subset with
[make_blast_subset.py](make_blast_subset.py)
`--truth <truth.tsv> --out-dir <dir> --strata <column> [<column> ...] --per-stratum <n> --seeds <seed> [<seed> ...]`;
it writes one ReadID list per seed for selecting both mates.

vsearch and minimap2 follow the same rules: same reference records, same
fastp-filtered FASTQs, settings frozen before examining accuracy, and the same
LCA conversion to the common taxonomy paths.
Query FASTAs must retain the `/1` and `/2` read ID suffixes so alignment hits can be assigned to mates.

## Evaluation

1. Use one conversion script, [convert_calls.py](convert_calls.py), that maps
   every method's calls into the same `ReadID`, `Taxonomy`, `status` schema.
   Its subcommands are `barcodemapper`, `best-hit`, `kraken2`, `blastn`,
   `vsearch` and `minimap2`, each taking `--read-ids`, `--input` and `--output`;
   `kraken2` also needs `--taxid-paths` and the alignment tools need
   `--accession-paths`. Its tests cover one-row-per-pair behavior,
   taxon-prefix mapping, and missing/unclassified cases on small fixtures.
2. On synthetic truth, report rank-specific precision, recall, wrong species
   overcalls, resolution depth, and false calls per million off-target pairs.
   Plot precision against recall across each method's *predeclared* settings;
   never optimize settings on the final confirmation set.
   [evaluate_calls.py](evaluate_calls.py) computes these:
   `--calls <method>=<calls.tsv>` (repeat per method) `--truth <truth.tsv>
   --out <dir>`, with an optional `--hashes-manifest` to verify that every
   method used the same inputs.
3. On real mocks, report taxa detected, unexpected calls, and fragment support
   under the same thresholds as Project 01. Keep reference-absent mock members
   visible in the denominator table. Use `--expected-taxa <expected-taxa.tsv>
   --run-accession <accession>` in place of `--truth`.
4. Report classifier-only wall time, CPU-hours, peak memory, and index size.
   Index construction time and memory are separate. Use identical CPU limits
   and the same compute-node class for timing comparisons; repeat on at least
   three independent job allocations to reflect cluster variation.
   [time_classifier.sh](time_classifier.sh) wraps one classifier command:
   `--label classifier-only|end-to-end --cpus <n> --input <fastq> --index <path>
   --output <timing.json> -- <command>`. Repeat `--input` and `--index` as needed.

## Deliverables and completion check

Save custom taxonomy conversion source and tests, FASTA record audit, Kraken2
index manifest, BLASTn configuration, converted calls, accuracy tables, timing
and resource logs, and all filled run records. The project is complete when
every method used the same FASTQ and reference record hashes for each dataset,
and the common evaluation script has verified those inputs.
