# BarcodeMapper Project 01 pilot: local evaluation

**Scope and caveats.** Preliminary 3-run pilot: one animal mock (COI, run ERR6001161), one fungal pooled-DNA sample (ITS, run SRR14351851) and one extraction blank (ERR6001173), each subsampled to 100k read pairs (200k reads) and run with the 'pilot-standard' BarcodeMapper settings. There is exactly one run per marker, so precision and recall are indicative only; no confidence intervals or between-replicate variation are available. Species recall is a lower bound: `synonym_status` is `unreviewed` for every expected taxon (exact paper binomial only, no synonyms invented), so a taxon reported under a current synonym is counted as a miss and its reported name as an unexpected call. The blank role (ERR6001173) is an inference in `sample-roles.tsv`.

## Method

- Counts in the taxonomy TSVs are **cumulative** (a taxon includes all its descendants); see `data_checks.tsv`. Exclusive counts are derived as parent minus sum of children.
- A taxon is *detected* at threshold t if its cumulative fragment count is >= t (t = 1, 2, 5). `Unclassified` is never a taxon.
- Matching: exact name, case-insensitive, `_` equals space. Expected genus = first word of the expected binomial; expected family = `source_families`.
- Precision = TP / detected at that rank (all detected taxa, including other kingdoms; a `kingdom_scope` row restricted to the target kingdom is also given).
- Recall has three denominators: `recall_all_expected` (all expected taxa at the rank), `recall_exact_species_db` (expected taxa derived from species with `coverage_status == exact_species`, i.e. the species is in the reference) and `recall_rank_db` (expected taxa with >0 reference records at that rank, from the genus/family record counts).

## Headline numbers (kingdom_scope = all)

| marker | rank | min frag | detected | expected (all / in-DB sp.) | TP | FP | precision | recall all | recall in-DB |
|---|---|---|---|---|---|---|---|---|---|
| COI | species | 1 | 8 | 28 / 26 | 5 | 3 | 0.62 | 0.18 | 0.19 |
| COI | species | 2 | 4 | 28 / 26 | 2 | 2 | 0.50 | 0.07 | 0.08 |
| COI | species | 5 | 0 | 28 / 26 | 0 | 0 | NA | 0.00 | 0.00 |
| COI | genus | 1 | 10 | 24 / 22 | 8 | 2 | 0.80 | 0.33 | 0.32 |
| COI | genus | 2 | 5 | 24 / 22 | 4 | 1 | 0.80 | 0.17 | 0.14 |
| COI | genus | 5 | 0 | 24 / 22 | 0 | 0 | NA | 0.00 | 0.00 |
| COI | family | 1 | 12 | 18 / 18 | 8 | 4 | 0.67 | 0.44 | 0.44 |
| COI | family | 2 | 6 | 18 / 18 | 3 | 3 | 0.50 | 0.17 | 0.17 |
| COI | family | 5 | 0 | 18 / 18 | 0 | 0 | NA | 0.00 | 0.00 |
| ITS | species | 1 | 12 | 38 / 33 | 5 | 7 | 0.42 | 0.13 | 0.15 |
| ITS | species | 2 | 8 | 38 / 33 | 2 | 6 | 0.25 | 0.05 | 0.06 |
| ITS | species | 5 | 3 | 38 / 33 | 1 | 2 | 0.33 | 0.03 | 0.03 |
| ITS | genus | 1 | 17 | 24 / 23 | 9 | 8 | 0.53 | 0.38 | 0.39 |
| ITS | genus | 2 | 8 | 24 / 23 | 3 | 5 | 0.38 | 0.12 | 0.13 |
| ITS | genus | 5 | 4 | 24 / 23 | 2 | 2 | 0.50 | 0.08 | 0.09 |
| ITS | family | 1 | 16 | 14 / 14 | 10 | 6 | 0.62 | 0.71 | 0.71 |
| ITS | family | 2 | 8 | 14 / 14 | 4 | 4 | 0.50 | 0.29 | 0.29 |
| ITS | family | 5 | 5 | 14 / 14 | 2 | 3 | 0.40 | 0.14 | 0.14 |

## Read fate

Pairs retained by fastp use `fastp_reads_after / 2` (assumes paired filtering). 'Matched' is `assigned_read_ids` from the manifest, which equals the sum of kingdom-level counts plus `Unclassified` (checked), i.e. matched fragments include those left `Unclassified`; that reading of `Unclassified` is inferred, not documented in the files.

| run | marker | pairs_input | pairs_retained_fastp | marker_matched_assigned_ids | frac_input_pairs_matched | resolved_species | resolved_genus | resolved_family | resolved_higher_order_to_kingdom | unclassified |
|---|---|---|---|---|---|---|---|---|---|---|
| ERR6001161 | COI | 100000 | 98155 | 21 | 0.00021 | 13 | 3 | 3 | 2 | 0 |
| ERR6001173 | none | 100000 | 98880 | 4 | 0.00004 | 0 | 0 | 1 | 3 | 0 |
| SRR14351851 | ITS | 100000 | 92383 | 226 | 0.00226 | 98 | 11 | 29 | 42 | 46 |

Fields available in the inputs: manifest `fastp_reads_before`, `fastp_reads_after`, `assigned_read_ids`; `*.results.json` holds only artifact paths/hashes (no counts). Per-fragment assignments and the pre-match-filter counts are not available locally.

## Surprises worth a look

- Very few pairs reach an accepted marker match (21, 226 and 4 of 100,000 pairs for the animal, fungal and blank runs); 100k pairs of shotgun data hold only a handful of marker fragments, so the species-level evidence per expected taxon is 0-3 fragments and precision/recall are dominated by sampling noise.
- The fungal run was mapped with the `all-fungi` reference yet received Animalia calls (26 fragments, mostly Cecidomyiidae, 21) and a Dipodidae call; check what non-fungal ITS records the all-fungi database contains.
- The blank (default reference) has Viridiplantae (Trebouxiophyceae, 1 fragment), Hymenoptera (2) and a Sminthurididae (1) fragment although it is a blank; Sminthurididae is also an expected family of the animal mock (*Sphaeridia pumilis*), so sample-to-blank carry-over cannot be excluded.
- The animal mock run received one Fungi/Polyporaceae fragment.

## Unexpected calls

- ERR6001161: 3 unexpected species, 6 fragments; 1 fragments assigned to a kingdom other than the target (Animalia). Full list: `unexpected_calls.tsv`.
- SRR14351851: 7 unexpected species, 87 fragments; 26 fragments assigned to a kingdom other than the target (Fungi). Full list: `unexpected_calls.tsv`.

Names that look like synonyms of expected taxa (flagged for manual review, **not** applied to any metric): `Scedosporium_prolificans` vs expected `Lomentospora prolificans`; `Nakaseomyces glabratus` vs expected `Candida glabrata`. `Eupelops_sp._SFH2024_01` is a placeholder species name in a genus whose expected species is `Eupelops occultus`. This is why species recall is a lower bound.

## Blank control (ERR6001173)

Kingdom-level calls: Animalia=3, Viridiplantae=1. Full lineage list in `blank_calls.tsv`. No precision/recall computed; calls are not interpreted as contamination without further evidence. Family-level calls also present in the animal-mock expected set are flagged in `in_animal_mock_expected_family`.

## Expected taxa without a derivable rank

- ERR6001161: genus not derivable for 0 taxa; family not derivable (empty `source_families`) for 2 taxa: Deuterosminthurus bicinctus, Poikilolaimus oxycerca. These are absent from the family-level expected set.
- SRR14351851: genus not derivable for 0 taxa; family not derivable (empty `source_families`) for 12 taxa: Candida albicans, Candida dubliniensis, Candida glabrata, Candida haemulonis, Candida parapsilosis, Candida tropicalis, Diutina catenulata, Diutina mesorugosa, Diutina rugosa, Geotrichum fermentans, Lomentospora prolificans, Yarrowia lipolytica. These are absent from the family-level expected set.

## Files

- `rank_metrics.tsv`
- `expected_detection.tsv`
- `unexpected_calls.tsv`
- `read_fate.tsv`
- `blank_calls.tsv`
- `data_checks.tsv`
- `plot_read_fate.svg`
- `plot_recall.svg`

matplotlib is not installed in this environment, so plots are SVG (this script has no matplotlib in this environment); the PNGs were rendered from them with `rsvg-convert -b white -z 2`. Install matplotlib and rerun to get native PNGs.
