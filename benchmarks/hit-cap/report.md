# Bowtie2 hit cap (`-k`) vs all hits (`-a`): 1M-pair benchmark

Run on Mjolnir on 2026-10-16/17 with BarcodeMapper at commit `9440fe0` (1.1.3 plus the opt-in `-k/--max-hits`), the `--standard` preset, 8 mapping threads. One replicate per condition.

## What was run

- **Inputs:** the first 1,000,000 read pairs of two deep public runs. Animal: ERR6001161 (soil invertebrate mock, COI) against the `default` reference. Fungal: SRR14351851 (pooled DNA of 44 fungal species, ITS) against the `all-fungi` reference. The 100k pilot inputs are prefixes of these.
- **Settings:** `-a` (all hits, the current default, used as the baseline) and `-k` 1, 5, 10, 50, 200. 12 pipeline runs in total.
- **Mapping-only probes:** Bowtie2 piped into `samtools view -b` on the animal reads, at 10k and 1M pairs, with and without `-a`.
- **Comparison:** per-read taxonomy of each capped run against the `-a` run (`benchmarks/mjolnir/compare_max_hits.py`). Result tables are `compare_animal.tsv`, `compare_fungal.tsv` and `resources.tsv` in this folder. Per-read tables are in `projects/05-hit-cap/work/<sample>/<setting>/taxonomy/` on Mjolnir.

## Cost: the cap removes almost all the mapping time

Mapping step only (the final successful attempt), 1M pairs:

| Sample | Setting | Mapping time | Peak memory |
|---|---|---|---|
| animal | `-a` | 3h05m (plus a failed 3 min attempt) | 52.8 GB |
| animal | `-k` 200 / 50 / 10 / 5 / 1 | 94 / 149 / 127 / 92 / 232 s | 22.0 / 22.1 / 23.2 / 21.8 / 23.2 GB |
| fungal | `-a` | 2h59m | 38.8 GB |
| fungal | `-k` 200 / 50 / 10 / 5 / 1 | 153 / 127 / 116 / 128 / 228 s | 23.5 / 23.4 / 23.3 / 23.3 / 23.3 GB |

- The whole pipeline took 10-14 min with a cap and 3h08m-3h17m with `-a`.
- With a cap, time is dominated by loading the index (about 1 min); caps from 5 to 200 are indistinguishable in cost. The ~22-23 GB is the index plus the one-time BAM header.

Mapping-only probe on the animal reads (no fastp, 8 threads):

| Pairs | Setting | Time | Bowtie2 / samtools peak | Alignment rows | BAM |
|---|---|---|---|---|---|
| 10k | no `-a` | 98 s | 18.5 / 6.5 GB | 20,000 | 213 MB |
| 10k | `-a` | 113 s | 19.2 / 6.5 GB | 68,867 | 213 MB |
| 1M | no `-a` | 88 s | 18.6 / 6.5 GB | 2,000,000 | 322 MB |
| 1M | `-a` | 8,144 s (2h16m) | 47.2 / 6.5 GB | 8,501,620 | 410 MB |

- Mapping 1M pairs without `-a` costs about the same as 10k (index loading dominates); `-a` is about 90 times slower at 1M.
- At 10k pairs `-a` is nearly free, so its cost is superlinear in read count, consistent with a minority of reads that match many records.

## Accuracy: effect of the cap on read assignments

Reads with an accepted hit under `-a`: 155 (animal), 2,335 (fungal). Counts are reads. Assignments are compared after reducing them to their resolved prefix (highest contiguous named ranks), which is how the pipeline's final table counts them. "Identical" means the same taxonomy as `-a`.

| Sample | Setting | Identical | Lost | More specific | Different lineage | Species-resolved (`-a`: 110 animal, 1,099 fungal) |
|---|---|---|---|---|---|---|
| animal | `-k` 200 | 154 (99.4%) | 1 | 0 | 0 | 109 |
| animal | `-k` 50 | 153 (98.7%) | 2 | 0 | 0 | 108 |
| animal | `-k` 10 | 141 (91.0%) | 4 | 10 | 0 | 111 |
| animal | `-k` 5 | 138 (89.0%) | 6 | 11 | 0 | 111 |
| animal | `-k` 1 | 108 (69.7%) | 18 | 29 | 0 | 118 |
| fungal | `-k` 200 | 2,331 (99.8%) | 0 | 4 | 0 | 1,100 |
| fungal | `-k` 50 | 2,324 (99.5%) | 1 | 10 | 0 | 1,099 |
| fungal | `-k` 10 | 2,251 (96.4%) | 2 | 82 | 0 | 1,128 |
| fungal | `-k` 5 | 2,121 (90.8%) | 22 | 192 | 0 | 1,184 |
| fungal | `-k` 1 | 1,423 (60.9%) | 182 | 726 | 4 | 1,566 |

- "Lost" = a read with an accepted hit under `-a` and none under `-k`; "more specific" = the capped run names a more specific taxon on the same lineage; "different lineage" = a lineage not on the `-a` path. No capped run was ever less specific than `-a`.
- Capping works by discarding hits, so the shared-taxonomy step sees less competing evidence and calls reads more specifically. At `-k` 1 that is 726 of 2,335 fungal reads, and species-level calls rise from 1,099 to 1,566. These extra calls come from less evidence, not from more information, so they are not better accuracy.
- Reads almost never move to a different lineage: 4 reads at `-k` 1 and none at any larger cap. The cost of a cap is mainly lost reads and over-specific calls.
- `-k` 50 and above leave 98.7-99.8% of reads unchanged and the species-level counts within 1 read.
- A few reads appear only with a cap (gained: 1-3 per fungal run). Bowtie2 `-k` uses a different search from `-a`, so the two are not strict subsets.

Accepted alignments per read under `-a` (the demand a cap has to meet):

| Sample | Reads | Alignments | Median | Max | Reads with >10 | >50 | >200 |
|---|---|---|---|---|---|---|---|
| animal | 155 | 1,782 | 2 | 248 | 29 | 11 | 1 |
| fungal | 2,335 | 24,057 | 2 | 3,020 | 302 | 49 | 14 |

## Accuracy against the expected taxa

The same scoring as the pilot (`benchmarks/pilot-01/evaluate_pilot.py`, expected taxa of the mock communities), applied to each setting; all rows, including other ranks, thresholds and kingdom scopes, are in `accuracy_by_setting.tsv`. A taxon counts as detected at >= 1 fragment; "recall" is over all expected taxa (COI: 28 species, 24 genera, 18 families; ITS: 38 species, 24 genera, 14 families). Species recall is a lower bound because synonyms were not reconciled (see `benchmarks/pilot-01/analysis/README.md`).

| Marker | Setting | Species TP / FP | Species precision | Species recall | Genus precision | Genus recall |
|---|---|---|---|---|---|---|
| COI (animal) | `-a` | 11 / 5 | 0.69 | 0.39 | 0.79 | 0.46 |
| COI (animal) | `-k` 200 / 50 | 11 / 5 | 0.69 | 0.39 | 0.79 | 0.46 |
| COI (animal) | `-k` 10 | 12 / 6 | 0.67 | 0.43 | 0.80 | 0.50 |
| COI (animal) | `-k` 5 | 13 / 5 | 0.72 | 0.46 | 0.80 | 0.50 |
| COI (animal) | `-k` 1 | 14 / 10 | 0.58 | 0.50 | 0.80 | 0.50 |
| ITS (fungal) | `-a` | 12 / 21 | 0.36 | 0.32 | 0.45 | 0.58 |
| ITS (fungal) | `-k` 200 / 50 | 12 / 21 | 0.36 | 0.32 | 0.45 | 0.58 |
| ITS (fungal) | `-k` 10 | 13 / 23 | 0.36 | 0.34 | 0.44 | 0.58 |
| ITS (fungal) | `-k` 5 | 14 / 25 | 0.36 | 0.37 | 0.42 | 0.58 |
| ITS (fungal) | `-k` 1 | 17 / 46 | 0.27 | 0.45 | 0.42 | 0.67 |

- **`-k` 50 and 200 give exactly the same true and false positives as `-a`** at species, genus and family level, for both markers. At the stricter threshold (>= 2 fragments) `-k` 50 differs from `-a` by one species on ITS.
- **Smaller caps recover a few more expected species but add false ones.** At `-k` 1 ITS gains 5 true species and 25 false ones (precision 0.36 to 0.27); COI gains 3 and 5 (0.69 to 0.58). Their extra recall is the over-specific calling described above, not more evidence.
- The absolute precision and recall are low for all settings, as they were in the pilot: the data hold few marker reads (a few hundred animal or a few thousand fungal fragments in 1M pairs), and the reference and synonym issues are unchanged from the pilot. Differences between settings are small counts, a handful of taxa each.

## Caveats

- **Single replicate, two samples.** "Identical to `-a`" measures agreement with the current default, not correctness; the expected-taxa scoring above checks correctness but on few taxa, with a synonym-limited recall.
- **The animal sample is thin:** 155 marker reads at 1M pairs, so its percentages move by about 0.6% per read.
- **`-k` does not choose the best N hits:** Bowtie2 stops after N valid alignments, which also count alignments later rejected by the `-m`/`-c` filters. That is why caps below 50 change calls for reads with many near-identical hits.
- **Resource formulas undersize `-a` runs:** the animal `-a` mapping job went out of memory at 46 GB and was retried at 92 GB; the fungal `-a` job finished 2 min 38 s inside its 3h02m time limit. Capped runs fit comfortably. The formulas scale with input size only, not with how many hits a read set produces.
- **Earlier species-count discrepancy, resolved:** the pilot analysis (98 fungal reads at species level) was right. A first version of the comparison script counted reads ending in empty-named ranks such as `g__;s__` as species-level (140). The pipeline keeps only the highest contiguous named ranks, so those reads count at family or above. The script now applies the same rule and reproduces the pilot's 98 / 11 / 29 (species / genus / family) exactly; the numbers above use the corrected rule.
- Times include cluster scheduling effects (different nodes); no timing replicates.

## Recommendation (from this one benchmark)

- Keep `-a` as the default for now and keep `-k` opt-in.
- A cap of 50-200 looks like a safe speed-up: about 99% or more of reads unchanged, the mapping step 70-120 times faster and peak memory 40-60% lower at 1M pairs. If a convenience preset is wanted, `--fast` = `-k 200` is the most conservative candidate.
- Do not offer caps of 10 or below as a normal setting: 3.6-9% of reads change assignment at `-k` 10 and 9-11% at `-k` 5, and `-k` 1 changes 30-39% of them, mostly by becoming over-specific.
- Before changing any default: (1) repeat on more samples and at 10M pairs, where the animal sample has enough marker reads; (2) try caps between 50 and 200 (for example 75, 100, 150), to see where, if anywhere, assignments start to differ from `-a`. Scoring against the expected taxa is done (above): caps of 50 and 200 matched `-a`.
