# Project 01: defined shotgun mock communities

## Question and scope

How often does BarcodeMapper detect known animal and fungal taxa in real
shotgun data, at what taxonomic rank, and how many unexpected taxa does it
report? This is a taxon-presence benchmark. A known member of a mock community
does not provide a ground-truth species label for every individual read.

## Candidate public datasets

| Marker | Dataset | Why use it | Database |
| --- | --- | --- | --- |
| COI | Soil invertebrate mixtures, ENA study [PRJEB45431](https://www.ebi.ac.uk/ena/browser/view/PRJEB45431) | Controlled mixtures of 28 soil invertebrate species, shotgun Illumina PE150 | `default.fa.gz` |
| ITS | Fungal pooled-DNA and pooled-biomass mixtures, SRA/ENA study [PRJNA725368](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA725368) | Two independently prepared mixtures with Illumina shotgun data; 43 expected Table 1 rows collapse to 38 distinct binomial names | `all-fungi.fa.gz` |

The animal [study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9170594/) and
fungal [study](https://journals.asm.org/doi/10.1128/mbio.02444-21) document
the member lists. For PRJNA725368, **select only the Illumina shotgun runs**;
the study also contains ITS amplicon and Nanopore runs. Save the paper's member
tables, dataset citations, and the exact run-selection rationale. Do not
assume every expected species has an eligible sequence in the frozen reference.

## Inputs to freeze

1. Obtain an ENA `read_run` file report for each study, saving the response as
   `inputs/ena-PRJEB45431.tsv` and `inputs/ena-PRJNA725368.tsv`. Include at
   least `run_accession`, `sample_title`, `library_strategy`,
   `library_layout`, `instrument_model`, `fastq_ftp`, `fastq_md5`, and
   `fastq_bytes`. The [ENA file-report API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)
   documents the accession query. Record the query URL and retrieval date.
2. Curate `inputs/selected-runs.tsv`: one row per included paired shotgun run,
   with study, run accession, sample type, reported member table, R1/R2 URLs,
   ENA MD5s, and inclusion decision. Check that the two FASTQ files are truly
   paired and correspond to shotgun sequencing. Keep excluded rows and reasons.
3. Download selected FASTQs to `inputs/raw/` and verify both ENA MD5s. Record
   SHA-256 values locally as well. Symlink them under `inputs/pairs/` with
   unique names ending `_1.fastq.gz` and `_2.fastq.gz`; do not merge biologically
   distinct preparations. Preserve accession-to-sample mapping.
4. Curate `inputs/expected-taxa.tsv` with sample, expected species, accepted
   synonym, marker, and whether the frozen database has an eligible species,
   genus, or only higher-rank reference. Decide taxonomy aliases before viewing
   BarcodeMapper outputs. Report reference-absent taxa separately.
   [annotate_expected_taxa.py](annotate_expected_taxa.py) makes the first exact-name
   coverage annotation after reference audits; review unresolved names and
   synonyms before computing recall.

## Runs

In the project's dedicated environment, run one sample at a time initially.
This keeps memory, scratch use, and failure diagnosis bounded. Use the same
preprocessing parameters and the standard mapping preset in the primary
analysis. Set the adapters to `auto` only if the libraries' adapters are not
known; record the choice. Example for an animal sample:

```bash
BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
PROJECT_ROOT="$BENCH_ROOT/projects/01-mock-communities"
barcodemapper --slurm \
  -1 "$PROJECT_ROOT/inputs/pairs/animal01_1.fastq.gz" \
  -2 "$PROJECT_ROOT/inputs/pairs/animal01_2.fastq.gz" \
  -d "$BENCH_ROOT/references/default.fa.gz" \
  -o "$PROJECT_ROOT/results/animal01-standard.tsv" \
  -t "$PROJECT_ROOT/work/animal01-standard" \
  --standard --adapter-r1 auto --adapter-r2 auto
```

For the fungal shotgun runs, switch to `all-fungi.fa.gz` and give every run a
separate result and work path. After the primary runs, rerun a predeclared
representative subset with `--ultrastrict`, `--strict`, and `--permissive` to
show the effect of stringency. Reuse identical input files and reference
snapshots; use a fresh work directory for every preset. Submit with the site's
approved Slurm controller pattern and retain controller and worker job IDs.

If extraction blanks or negative libraries from the same study are available,
include them. If not, state that biological false-positive rates in extraction
blanks cannot be estimated from these public datasets; use the known off-target
simulations in Project 02 as a separate specificity test.

## Evaluation

- Primary: taxon-level precision and recall at species, genus, and family;
  count a taxon as detected at thresholds of at least 1, 2, and 5 assigned
  fragments. Use expected taxa in the database and all expected taxa as two
  clearly labelled recall denominators.
- Report unexpected taxon names and fragment support, including whether they
  are plausible contaminants or conflicts with published mock composition.
  A taxon outside the supplied member list is an *unexpected call*, not proof
  of a laboratory false positive.
- Report the fraction of input pairs retained by fastp, pairs with an accepted
  marker match, and matched pairs resolved to species/genus/higher ranks or
  `Unclassified`. Plot COI and ITS separately.
- Secondary: association between input DNA or biomass and marker-matching
  fragment counts, with no organism-abundance claim.

## Deliverables and completion check

Save accession metadata, verified FASTQs or an immutable local location,
member tables, per-run manifests, per-fragment assignments, a rank-by-rank
results table, unexpected-call audit, plots, and filled run records. The
project is complete when every selected run has a verified checksum, a
successful BarcodeMapper manifest, and a traceable expected-taxon entry.
