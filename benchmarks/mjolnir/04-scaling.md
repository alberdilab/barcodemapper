# Project 04: runtime and storage scaling

## Question

How do BarcodeMapper wall time, CPU use, memory, and temporary disk grow with
read count, marker-hit rate, and sample count on Mjolnir? Bowtie2's all-hit
mode can produce many alignments per fragment, so temporary BAM size is a
first-class outcome, not merely a nuisance. The [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
notes that `-a` can be very slow for repetitive references.

## Fixed inputs and planned matrix

Use the frozen `default.fa.gz` reference for the primary curve. Create one
background-heavy real shotgun stream from verified public reads, and one
marker-rich stream from Project 02's synthetic target reads whose taxa are
represented in `default.fa.gz`. Unique read
IDs are required across concatenated sources. Save scripts, selected accessions,
order, seeds, pair counts, and SHA-256 hashes for every derived FASTQ.

| Axis | Values | Purpose |
| --- | --- | --- |
| Input pairs | Pilot 0.1M, then 1M, 10M, 50M if feasible | Read-count growth |
| Marker-hit fraction | Background-heavy and marker-rich | Alignment multiplicity |
| Sample layout | One sample and ten samples with the same total pair count | Scheduler and per-sample overhead |
| Repetition | Three runs on distinct allocations per condition | Run-to-run variability |

Build nested input sizes from one deterministic ordering; each smaller input
must be a prefix or subset of the larger input. Do not achieve the 50M-pair
condition by repeatedly copying the same small read set. If 50M independent
pairs cannot be assembled, reduce the ceiling and document it. Keep read
length, fastp settings, reference, and BarcodeMapper preset fixed at the
primary `--standard` condition. A secondary curve may use `--strict` after the
primary results are frozen.

## Pilot and resource gate

1. Run the 0.1M pilot first with a new `-t` directory. Record input bytes,
   pairs retained by fastp, accepted marker-matching pairs, alignment rows,
   Bowtie2 BAM bytes, peak scratch bytes, stage times, and Slurm resources.
2. Estimate scratch demand for 1M, 10M, and 50M using **observed bytes per
   input pair** for both hit-rate conditions, with room for intermediate files
   and retries. Confirm Mjolnir quotas and available partition limits before
   submitting larger jobs. A pilot failure or non-linear alignment explosion
   is itself a result to report, not a reason to silently drop that stratum.
3. Freeze the job resource request and Slurm profile for the three repeats of
   each condition. The packaged profile caps threads at eight; record the
   effective threads, job concurrency, memory limits, and any site-specific
   profile changes. Do not tune resources separately after seeing individual
   replicate timings unless the entire condition is rerun.

## Runs and measurement

Use the same pattern as Project 01, with one work and output path per
`condition/replicate`. Use the prebuilt reference index for analysis so index
creation is not folded into the input-size curve. Measure and publish the
database build/index cost separately from a cold build. Keep all intermediate
files until their sizes and checksums are recorded.

For each run record:

- Slurm submit, start, and end time; queue wait and execution wall time
  separately. Aggregate worker-job CPU-hours and MaxRSS from `sacct`; the
  Snakemake controller alone does not measure the work it submits.
- Per-stage elapsed time for fastp, Bowtie2 plus samtools collate, extraction,
  taxonomy assignment, aggregation, and HTML reporting. Capture job IDs and
  stage logs from the frozen Snakemake run.
- Peak scratch occupancy, input/output bytes, BAM bytes, extracted TSV bytes,
  final report bytes, and index bytes. Sample scratch occupancy at a fixed
  interval or use filesystem accounting; an end-of-run `du` is not a peak.
- Input pairs/s and accepted marker-matching pairs/s, plus alignment records
  per marker-matching pair. Record node model and filesystem location.

The one-sample/ten-sample comparison must keep **total** pairs constant.
Report maximum concurrency and total allocated core-hours, since ten samples
can finish sooner while consuming more cluster resources. Prefer paired
conditions on the same node class; do not subtract queue time from the raw
timestamp without recording it separately.

## Analysis and completion check

Plot wall time, CPU-hours, MaxRSS, peak scratch, and BAM size against input
pairs for both hit-rate strata. Show all repeats, not only the median, and
identify any incomplete or resource-limited runs. Report the largest
successfully completed dataset rather than extrapolating beyond it. The
project is complete when all attempted conditions have a run record,
successful outputs have checksum manifests, and failed conditions retain
their Slurm reason and last complete stage.

## Pilot result: background-heavy, 100k pairs, rep1 (clean-source run)

Controller job 45903725, BarcodeMapper 1.1.2, `default.fa.gz`
(sha256 `462c2cc6…3887`), `--standard`, run 2026-10-13 on Mjolnir.

| Quantity | Observed |
| --- | --- |
| Input pairs / reads after fastp | 100,000 / 196,310 reads |
| Accepted marker-matching pairs | 21 |
| Controller wall time | 2 h 07 min (1 CPU, 4 GB) |
| Bowtie2 mapping job (45903732) | 1 h 48 min, 8 CPUs, 3 h 53 min CPU |
| Bowtie2 mapping MaxRSS | about 199 GB (304 GB requested) |
| Downstream jobs (extract, aggregate, report) | 4-11 s each, under 1 MB RSS |
| Work directory after run | 17 GB |

Caveats:

- The input is the ERR6001161 animal pilot pairs, so `rep1.tsv` is byte-identical
  to `ERR6001161-pilot-standard.tsv`. It is not an independent background
  stream, and there is no replicate yet.
- `rep1.scratch.csv` is empty: the scratch sampler job (45903726) failed
  because `sample_scratch.py` used `capture_output`, unsupported by the
  system Python 3.6. Fixed by using `stdout/stderr=PIPE` and
  `universal_newlines`; peak scratch is therefore not measured for this run.
- Mapping memory of about 199 GB for 100k pairs is the open resource-gate
  item. Check whether all-hit mode against the full reference is the cause
  before extrapolating to 1M, 10M or 50M pairs.
