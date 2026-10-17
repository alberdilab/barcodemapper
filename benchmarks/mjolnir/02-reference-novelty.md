# Project 02: reference novelty and barcode ambiguity

## Question

When an exact source taxon is missing from the reference, does BarcodeMapper
return the deepest defensible ancestor instead of overcalling a related
species? How do identical barcodes with conflicting taxonomy affect the
resolution of paired shotgun reads? This project supplies read-pair truth,
which Project 01 cannot supply for every real read.

## Frozen inputs and design

Use the retained BOLD and UNITE source FASTAs from the [shared setup](README.md)
and the `all-fungi.fa.gz` taxon scope. Create separate COI and ITS test sets;
do not pool their metrics. Freeze a `design.tsv` containing every selected
source accession, marker, complete source taxonomy, sequence hash, sequence
length, and stratum. Selection must be deterministic from a recorded seed and
must happen **before** inspecting model outputs.

Define four disjoint strata for each marker, with several families and genera
represented where the source data allow:

1. **Present species:** source species and independent conspecific records
   remain in the training reference. Do not simulate from the exact FASTA
   record still indexed; use a held-out accession where possible.
2. **Absent species:** remove *all* source-species records from the training
   reference, retaining other species from its genus.
3. **Absent genus:** remove all records from the source genus, retaining
   representatives of the same family.
4. **Ambiguous barcode:** choose identical or near-identical sequence groups
   that carry conflicting species or genus labels; preserve those conflicts in
   the training reference and list all involved taxonomies.

Group selected source taxa into a small, predeclared number of folds (for
example, five), with one indexed reference per fold and condition. This avoids
rebuilding a multi-gigabyte index for every taxon. In the absent-genus folds,
verify that each held-out genus still has a represented family in the index.

Also generate a negative/background stratum from independently sourced
non-target shotgun reads or genomes. Record their accessions and hashes. Do
not use uniformly random DNA as the only negative control; it provides a weak
test of off-target mapping. Screen inputs for adapters and duplicate read IDs.

## Read generation and reference variants

- Partition at the **taxon** level for species/genus holdouts. Verify no
  accession from a held-out taxon remains in its training FASTA. Audit exact
  sequence matches that cross taxon labels and report them separately.
- Construct immutable marker-specific FASTAs for each fold and condition from the frozen
  source records. Maintain BarcodeMapper's `accession|marker|taxonomy` header
  format. Build matching Bowtie2 indexes once and retain the FASTA, index,
  partition table, and their hashes.
- Simulate paired reads with pinned ART 2016.06.05, a documented
  fragment-length distribution and its `HS25` quality profile (primary PE150;
  optional PE75 for short ITS records), and a deterministic seed. ART's
  [documentation](https://www.niehs.nih.gov/research/resources/software/biostatistics/art)
  gives the paired-end options. For example, the PE150 condition can start
  with `art_illumina -ss HS25 -i heldout.fa -p -l 150 -f 20 -m 350 -s 30 -rs 1337 -sam -o reads_`;
  adjust the insert distribution only before the design is frozen and record
  it. Use ART's SAM truth to write one row per
  fragment: read ID, source accession, marker, source taxonomy, stratum,
  fragment coordinates, and seed. Do not generate PE150 fragments from
  sequences too short to support the specified insert length; report the
  number excluded by marker and taxon.
- Mix known target and off-target pairs at predeclared proportions, including
  low target fractions representative of shotgun data. Use unique read IDs
  across all sources. Save the exact FASTQs and SHA-256 hashes.

The exact ART command, profile, and version become part of the run record
**before** data generation. The source of any helper for paired FASTQ mixing
and truth-table construction must be committed and its SHA-256 recorded.

## Runs

Run each test FASTQ against its matching immutable reference variant with
`--standard` first. Then run all four presets on the same inputs and database
to obtain the stringency tradeoff. Use a separate `-t` working directory for
each marker, holdout condition, seed, and preset. Do not select a preset on
the final test set; reserve at least one seed or taxon set as a locked
confirmation set.

Example naming convention:

```text
inputs/PE150/COI/species-holdout/seed-1337_1.fastq.gz
inputs/PE150/COI/species-holdout/seed-1337_2.fastq.gz
references/COI/species-holdout.fa.gz
results/COI/species-holdout/seed-1337/standard.tsv
work/COI/species-holdout/seed-1337/standard/
```

The BarcodeMapper invocation uses the same `-1`, `-2`, `-d`, `-o`, `-t`, and
`--slurm` structure as [Project 01](01-mock-communities.md). Its raw per-pair
assignments are in `work/.../taxonomy/<sample>.tsv`.

## Evaluation

- Join truth and predictions by **read ID**, not by taxonomy table row.
  Independently verify one output row at most per pair.
- For each taxonomic rank, report precision, recall, and fraction of pairs
  resolved at that rank. A correct ancestor is a valid higher-rank assignment
  and a lower-rank abstention. An assignment outside the true lineage is a
  wrong call; an unsupported species call is explicitly an overcall.
- For absent-species and absent-genus strata, plot the deepest correct rank
  and wrong-overcall frequency against sequence distance to the nearest
  retained reference. Report exact-duplicate cross-label cases separately.
- For negative/background pairs, report accepted marker hits and false taxon
  calls per million input pairs. Include unclassified and no-hit fractions.
- Bootstrap uncertainty by source taxon or simulated sample/seed, not by
  treating millions of fragments from one source as independent taxa.

## Deliverables and completion check

Archive `design.tsv`, frozen source/negative accessions and hashes, partition
audit, simulator source or version, truth table, FASTQ hashes, reference
variant hashes, per-run manifests, per-pair predictions, and analysis plots.
The project is complete only after an automated audit confirms no source
taxon intended for holdout appears in its corresponding training FASTA.
