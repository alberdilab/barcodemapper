# Bowtie2 hit cap (`-k`): extension to more samples and intermediate caps

Run on Mjolnir on 2026-10-17 with BarcodeMapper at commit `9440fe0` (1.1.3 plus the opt-in `-k/--max-hits`), the `--standard` preset, 8 mapping threads, the first 1,000,000 read pairs of each run. One replicate per condition. This extends the first benchmark (`report.md`: caps 1, 5, 10, 50, 200 against `-a` on one animal and one fungal sample) with the two follow-ups it recommended.

## What was run

- **Intermediate caps:** `-k` 75, 100 and 150 on the existing animal (ERR6001161, COI, `default` reference) and fungal (SRR14351851, ITS, `all-fungi` reference) inputs.
- **More samples, `-a` baseline plus `-k` 50 and 200:** three more animal mock communities (ERR6001164 = MC2_1, ERR6001167 = MC3_1, ERR6001170 = MC4_1) and one more fungal sample (SRR14351849, pooled biomass). That is 12 pipeline runs, plus 6 for the caps, 18 in total. All 19 jobs (including the input-cutting job) finished with exit code 0.
- **Analysis:** `analyze_extension.py` in this folder reads `data/` (result tables and per-read taxonomy pulled from Mjolnir) and writes `extension_compare.tsv` (per-read comparison with `-a`, same rules as `report.md`) and `extension_accuracy.tsv` (scoring against the expected taxa, same code as the pilot). Run times and memory come from `sacct` for the Bowtie2 mapping job of each run.

## Cost: the cap again removes almost all mapping time

Mapping step only (final successful attempt), 1M pairs:

| Sample | Setting | Mapping time | Peak memory |
|---|---|---|---|
| ERR6001164 (animal) | `-a` | 37.7 min | 54.8 GB |
| ERR6001164 | `-k` 50 / 200 | 114 / 213 s | 23.8 / 24.0 GB |
| ERR6001167 (animal) | `-a` | 2h45m | 55.2 GB |
| ERR6001167 | `-k` 50 / 200 | 112 / 101 s | 23.8 / 22.7 GB |
| ERR6001170 (animal) | `-a` | 5h13m | 53.5 GB |
| ERR6001170 | `-k` 50 / 200 | 108 / 122 s | 22.6 / 23.9 GB |
| SRR14351849 (fungal) | `-a` | 2h18m | 42.7 GB |
| SRR14351849 | `-k` 50 / 200 | 101 / 131 s | 22.6 / 24.2 GB |
| animal (ERR6001161) | `-k` 75 / 100 / 150 | 113 / 204 / 114 s | 23.8 / 23.9 / 23.9 GB |
| fungal (SRR14351851) | `-k` 75 / 100 / 150 | 134 / 133 / 135 s | 24.0 / 24.0 / 24.1 GB |

- Across the four new samples, capped runs take 1.7-3.6 min and `-a` takes 38 min to 5 h 13 min: 20 to 170 times slower. The spread among the `-a` runs (38 min to 5 h for similar read counts) shows that `-a` cost depends on how many repetitive reads a sample has, which cannot be predicted from its size.
- With a cap, time is the index load (about 1-2 min) and does not depend on the cap value, as in the first benchmark. Peak memory is 22-24 GB at every cap, against 43-55 GB with `-a`.

### Every `-a` run ran out of memory on its first attempt

| Sample | First attempt | Retry |
|---|---|---|
| ERR6001164 | out of memory at 46 GB (3 min) | 92 GB: finished in 38 min |
| ERR6001167 | out of memory at 46 GB (3.5 min) | 92 GB: **timed out at the 5h58m limit**; then 184 GB with 11h56m: finished in 2h45m |
| ERR6001170 | out of memory at 46 GB (3.5 min) | 92 GB: finished in 5h13m |
| SRR14351849 | out of memory at 47 GB (56 min) | 93 GB: finished in 2h18m |

- The measured peak is 43-55 GB, so the retry at about twice the request was enough in 3 of 4 cases. For ERR6001167 the 92 GB retry was not the problem: it was killed by the time limit, which the pipeline's resource formulas set from input size only. The whole `-a` pipeline for ERR6001167 took 11h59m including these failures.
- This confirms the caveat in `report.md`: resource formulas undersize `-a` runs. Capped runs fit in the first request every time (all 24 capped mapping jobs across both reports, none retried).

## Accuracy: effect of the cap on read assignments

Counts are reads whose resolved taxonomy equals the `-a` run ("identical"); "lost" has an accepted hit under `-a` and none under the cap; "more specific" is a cap-only refinement on the same lineage. No capped run was ever less specific than `-a`, and no read moved to a different lineage at any cap of 10 or above.

### New samples, `-k` 50 and 200

| Sample | Reads under `-a` | `-k` | Identical | Lost | More specific | Gained | Species-resolved (`-a` → cap) |
|---|---|---|---|---|---|---|---|
| ERR6001164 | 187 | 50 | 184 (98.4%) | 3 | 0 | 0 | 127 → 125 |
| ERR6001164 | 187 | 200 | 185 (98.9%) | 1 | 1 | 0 | 127 → 126 |
| ERR6001167 | 97 | 50 | 96 (99.0%) | 0 | 1 | 1 | 53 → 54 |
| ERR6001167 | 97 | 200 | 97 (100%) | 0 | 0 | 1 | 53 → 54 |
| ERR6001170 | 111 | 50 | 111 (100%) | 0 | 0 | 0 | 72 → 72 |
| ERR6001170 | 111 | 200 | 111 (100%) | 0 | 0 | 0 | 72 → 72 |
| SRR14351849 | 2,720 | 50 | 2,713 (99.7%) | 1 | 6 | 0 | 1,425 → 1,426 |
| SRR14351849 | 2,720 | 200 | 2,719 (99.96%) | 0 | 1 | 0 | 1,425 → 1,425 |

- `-k` 200 leaves 98.9-100% of reads unchanged on every new sample; `-k` 50 leaves 98.4-100%.
- The animal samples are thin (97-187 marker reads), so one read is 0.5-1%. The only sample with enough reads to give a stable percentage is the fungal one (2,720 reads): 99.7% at `-k` 50, 99.96% at `-k` 200.

### Intermediate caps on the first two samples

| Sample | `-k` | Identical | Lost | More specific | Gained | Species-resolved (`-a`: 110 animal, 1,099 fungal) |
|---|---|---|---|---|---|---|
| animal | 50 | 153 (98.7%) | 2 | 0 | 0 | 108 |
| animal | 75 | 153 (98.7%) | 2 | 0 | 0 | 108 |
| animal | 100 | 153 (98.7%) | 2 | 0 | 0 | 108 |
| animal | 150 | 154 (99.4%) | 1 | 0 | 0 | 109 |
| animal | 200 | 154 (99.4%) | 1 | 0 | 0 | 109 |
| fungal | 50 | 2,324 (99.5%) | 1 | 10 | 1 | 1,099 |
| fungal | 75 | 2,330 (99.8%) | 1 | 4 | 1 | 1,098 |
| fungal | 100 | 2,330 (99.8%) | 1 | 4 | 1 | 1,098 |
| fungal | 150 | 2,331 (99.8%) | 0 | 4 | 1 | 1,099 |
| fungal | 200 | 2,331 (99.8%) | 0 | 4 | 1 | 1,100 |

- Between 50 and 200 the change is gradual and small. On the fungal sample the step from 50 to 75 removes 6 of the 10 over-specific reads and nothing changes from 75 to 200 except one lost read recovered at 150. On the animal sample, 2 reads are lost at 50-100 and 1 at 150-200.
- The 1-2 reads lost on the animal sample at every cap up to 200 are not recovered by raising the cap within this range; I did not look at them individually; `report.md` notes that Bowtie2 `-k` searches differently from `-a`.

## Accuracy against the expected taxa

Same scoring as the pilot and `report.md` (detected at >= 1 fragment; recall over all expected taxa; species recall is a lower bound because synonyms were not reconciled). All rows, including other ranks, thresholds and kingdom scopes, are in `extension_accuracy.tsv`.

| Sample | Setting | Species TP / FP | Species precision | Species recall | Genus TP / FP |
|---|---|---|---|---|---|
| ERR6001164 | `-a`, `-k` 50, `-k` 200 | 8 / 4 | 0.67 | 0.29 | 11 / 3 |
| ERR6001167 | `-a`, `-k` 50, `-k` 200 | 8 / 3 | 0.73 | 0.29 | 10 / 2 |
| ERR6001170 | `-a`, `-k` 50, `-k` 200 | 7 / 3 | 0.70 | 0.25 | 10 / 2 |
| SRR14351849 | `-a`, `-k` 50, `-k` 200 | 12 / 22 | 0.35 | 0.32 | 13 / 19 |
| animal | `-a`, `-k` 50, 75, 100, 150, 200 | 11 / 5 | 0.69 | 0.39 | 11 / 3 |
| fungal | `-a`, `-k` 50, 75, 100, 150, 200 | 12 / 21 | 0.36 | 0.32 | 14 / 17 |

- **At >= 1 fragment, `-k` 50 to 200 give exactly the same true and false positives as `-a`** at species, genus and family level, on all six samples.
- At stricter thresholds a few taxa change. `-k` 150 and 200 differ from `-a` in no row on any sample. `-k` 50 differs in 8 rows: ERR6001164 loses 1 species and 1 genus at >= 5 fragments; the fungal sample and SRR14351849 gain one species at >= 2 fragments (and the fungal sample loses a false one). `-k` 75 and 100 differ in 2 rows on the fungal sample only (one false positive less at >= 2 fragments).
- The absolute precision and recall are low for every setting, as in the pilot (few marker reads, reference gaps, unreconciled synonyms), so none of these differences should be read as better or worse accuracy.

## Separate micro-benchmark: skipping the BAM step

`benchmarks/bam-skip/` compares the post-Bowtie2 work on the existing `-a` BAMs: `samtools view -b` plus pysam (current), a streaming SAM extractor (`sam_py`), pysam on SAM text, and an awk floor. Results are in `data/results/bam-skip/`.

- **Fungal** (2.64M alignments, 3 replicates): the current path takes 52.8-54.4 s and 6.5 GB; `sam_py` takes 5.2-5.3 s and 13 MB; pysam on SAM 46.9-50.1 s and 6.4 GB; awk 13.9-15.0 s and 4.5 MB. All four variants reproduce the current output exactly (24,057 lines, `equivalence.tsv`).
- **Animal:** `timings.tsv` has only a header and no timing rows; only `sizes.tsv` was written (2.14M alignments, 2.0 GB SAM, 323 MB BAM). So there is no animal timing result and no equivalence check for it; I did not investigate why.
- These timings leave out Bowtie2 and BAM writing itself, so the saving for a full run is on the order of a minute, which matters only for capped runs.

## Caveats

- **Single replicate per condition.** "Identical to `-a`" measures agreement with the current default, not correctness. The expected-taxa scoring checks correctness but covers few taxa.
- **Animal samples have 97-187 marker reads at 1M pairs.** Their percentages move by about 0.5-1% per read, so the conclusion for animals rests on the number of lost reads (0-3), not on the percentages.
- **Per-read hit counts were not collected for the new samples**, so I cannot say how many reads exceed 50 or 200 hits there. The first benchmark found 11 of 155 (animal) and 49 of 2,335 (fungal) reads above 50.
- **Peak memory for `-a` retries is the measured RSS,** not the request. The 184 GB request for ERR6001167 was a retry multiplier, not the need.
- Times include different compute nodes; there are no timing replicates.

## Recommendation (updated from both benchmarks)

- Keep `-a` as the default for now and keep `-k` opt-in.
- **`-k` 200 is now the best-supported candidate for a `--fast` preset.** Across six samples it matches the `-a` scoring against the expected taxa in every row, changes 0-1.1% of reads on the animal samples (0.5-1% per read, 1 read at most) and 0.04-0.2% on the two fungal ones, and cuts the mapping step from 38 min-5 h to under 4 min and peak memory by about half. `-k` 150 is equivalent on the two samples it was run on.
- `-k` 50-100 also looks safe in cost and read agreement (98.4-99.8%) but is the first setting where the stricter-threshold scoring starts to move; there is no cost reason to prefer it over 200, since time does not depend on the cap above about 5.
- Do not offer caps of 10 or below as a normal setting (see `report.md`).
- **Before changing a default:** repeat at 10M pairs on a sample with more marker reads, and count reads with many hits on the new samples. Separately, fix the resource formulas for `-a`: all four first attempts ran out of memory and one retry also hit the time limit.
