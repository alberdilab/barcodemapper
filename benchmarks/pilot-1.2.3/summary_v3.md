# Pilots 1.2.3: comparator evaluation on real mocks, v3 (exploratory)

Methods: BarcodeMapper 1.2.3 `--standard` (-k 200), `--all-hits`, single best-hit ablation, vsearch (n_mismatch rerun), minimap2 preset (primary minimap2 row), minimap2 sensitive (`-f 20000,100000`, labelled variant), Kraken2 (confidence 0.0, primary) and Kraken2 (confidence 0.1, predeclared variant). The first vsearch run is invalid (N-run artifact) and appears only in the appendix. BLASTn is still pending. Supersedes summary.md (v1) and v2.

Truth: Project 01 `expected-taxa.tsv` (exact paper binomials). Fragments = read pairs. Table cells: detected/expected taxa at min 1 fragment (min 2 fragments in parentheses). Unexpected = called taxa not expected, with total fragments in parentheses. `classified fraction` = classified pairs / retained pairs. `pairs with hit`: pairs not `no_hit` (for Kraken2, which has no no_hit status, classified pairs). FASTQ and reference sha256 were verified for every method in each run (identical across methods).

## SRR14351851: fungal pooled DNA, all-fungi reference

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 227 | 181 | 0.00196 | 5/38 (2) | 9/24 (4) | 10/14 (4) | 8 (88) | 8 | 6 |
| barcodemapper_allhits | 226 | 180 | 0.00195 | 5/38 (2) | 9/24 (3) | 10/14 (4) | 7 (87) | 8 | 6 |
| best_hit | 227 | 227 | 0.00246 | 11/38 (7) | 13/24 (9) | 12/14 (8) | 17 (115) | 9 | 7 |
| vsearch (n_mismatch) | 446 | 401 | 0.00434 | 13/38 (10) | 15/24 (13) | 13/14 (12) | 25 (170) | 26 | 22 |
| minimap2_preset | 712 | 577 | 0.00625 | 3/38 (2) | 4/24 (2) | 6/14 (4) | 18 (186) | 25 | 23 |
| minimap2_sensitive | 716 | 569 | 0.00616 | 3/38 (2) | 4/24 (2) | 6/14 (4) | 18 (185) | 25 | 23 |
| Kraken2 (conf 0.0) | 7092 | 7092 | 0.07677 | 18/38 (12) | 16/24 (12) | 13/14 (13) | 1312 (2014) | 1652 | 1201 |
| Kraken2 (conf 0.1) | 1061 | 1061 | 0.01149 | 7/38 (5) | 10/24 (7) | 12/14 (8) | 13 (289) | 15 | 12 |

Top unexpected species (fragments):

- barcodemapper_standard: Zygotorulaspora_dagestanica (58); Cordyceps_fumosorosea (17); Nakaseomyces_glabratus (3); Scedosporium_prolificans (3)
- barcodemapper_allhits: Zygotorulaspora_dagestanica (58); Cordyceps_fumosorosea (17); Nakaseomyces_glabratus (3); Scedosporium_prolificans (3)
- best_hit: Zygotorulaspora_dagestanica (59); Cordyceps_fumosorosea (18); Vishniacozyma_pseudopenaeus (8); Holtermannia_saccardoi (4)
- vsearch (n_mismatch): Zygotorulaspora_dagestanica (78); Cordyceps_fumosorosea (20); Tuber_oregonense (17); Nakaseomyces_glabratus (8)
- minimap2_preset: Zygotorulaspora_dagestanica (125); Cordyceps_fumosorosea (21); Chrysomycena_dunicola (7); Holtermannia_saccardoi (7)
- minimap2_sensitive: Zygotorulaspora_dagestanica (125); Cordyceps_fumosorosea (20); Chrysomycena_dunicola (7); Holtermannia_saccardoi (7)
- Kraken2 (conf 0.0): Zygotorulaspora_dagestanica (330); Cordyceps_fumosorosea (50); Chrysomycena_dunicola (41); Vishniacozyma_pseudopenaeus (26)
- Kraken2 (conf 0.1): Zygotorulaspora_dagestanica (224); Cordyceps_fumosorosea (40); Vishniacozyma_pseudopenaeus (7); Tulostoma_eckbladii (5)

## ERR6001161: animal mock, default reference

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 21 | 21 | 0.00021 | 5/28 (2) | 8/24 (4) | 8/18 (3) | 3 (6) | 2 | 4 |
| barcodemapper_allhits | 21 | 21 | 0.00021 | 5/28 (2) | 8/24 (4) | 8/18 (3) | 3 (6) | 2 | 4 |
| best_hit | 21 | 21 | 0.00021 | 5/28 (2) | 9/24 (4) | 9/18 (3) | 6 (10) | 2 | 4 |
| vsearch (n_mismatch) | 24 | 24 | 0.00024 | 5/28 (2) | 8/24 (5) | 8/18 (4) | 5 (9) | 2 | 4 |
| minimap2_preset | 78 | 74 | 0.00075 | 4/28 (2) | 7/24 (4) | 8/18 (4) | 4 (20) | 4 | 7 |
| minimap2_sensitive | 79 | 75 | 0.00076 | 4/28 (2) | 7/24 (4) | 8/18 (4) | 4 (20) | 4 | 8 |
| Kraken2 (conf 0.0) | 724 | 724 | 0.00738 | 5/28 (2) | 8/24 (5) | 8/18 (5) | 169 (216) | 219 | 243 |
| Kraken2 (conf 0.1) | 74 | 74 | 0.00075 | 5/28 (2) | 8/24 (5) | 8/18 (4) | 3 (5) | 3 | 7 |

Top unexpected species (fragments):

- barcodemapper_standard: Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Steganacarus_spinosus (1)
- barcodemapper_allhits: Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Steganacarus_spinosus (1)
- best_hit: Eupelops_sp._SFH2024_01 (3); Poikilolaimus_oxycercus (2); Seinura_italiensis (2); Enchytraeus_bigeminus_L1 (1)
- vsearch (n_mismatch): Eupelops_sp._SFH2024_01 (3); Poikilolaimus_oxycercus (2); Seinura_italiensis (2); Enchytraeus_dudichi (1)
- minimap2_preset: Hemiberlesia_rapax (14); Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Enchytraeus_sp._CEF (1)
- minimap2_sensitive: Hemiberlesia_rapax (14); Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Enchytraeus_sp._CEF (1)
- Kraken2 (conf 0.0): Stereopsis_radicans (10); Chlamydomonas_reinhardtii (8); Moniliophthora_canescens (4); Eupelops_sp._SFH2024_01 (3)
- Kraken2 (conf 0.1): Eupelops_sp._SFH2024_01 (2); Stereopsis_radicans (2); Seinura_italiensis (1)

## ERR6001173: extraction blank (negative control), default reference

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |
| barcodemapper_allhits | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |
| best_hit | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 2 (2) | 2 | 2 |
| vsearch (n_mismatch) | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 1 (1) | 1 | 2 |
| minimap2_preset | 20 | 18 | 0.00018 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 1 (1) | 1 | 2 |
| minimap2_sensitive | 21 | 19 | 0.00019 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |
| Kraken2 (conf 0.0) | 1140 | 1140 | 0.01153 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 247 (287) | 345 | 367 |
| Kraken2 (conf 0.1) | 19 | 19 | 0.00019 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 2 (2) | 2 | 4 |

Top unexpected species (fragments):

- barcodemapper_standard: none
- barcodemapper_allhits: none
- best_hit: Deuterostichococcus_epilithicus (1); Sphaeridia_pumilis (1)
- vsearch (n_mismatch): Deuterostichococcus_epilithicus (1)
- minimap2_preset: Chironomus_riparius (1)
- minimap2_sensitive: none
- Kraken2 (conf 0.0): Carinina_chocolata (4); Acartia_hudsonica (3); Buffonellaria_biavicularis (3); Chlamydomonas_reinhardtii (3)
- Kraken2 (conf 0.1): Deuterostichococcus_epilithicus (1); Stereopsis_radicans (1)

## Exploratory timings (not the final timing runs)

Classifier-only wall time, CPU-hours and peak RSS from `time_classifier.sh` (`timing.json`), 16 CPUs, on mixed node classes (compute nodes comp01 to comp16 and GPU node gpu03 of the same partition, so wall times are not strictly comparable). Index or database load is included. vsearch and minimap2 ran R1 and R2 separately; Kraken2 ran the pair in one command.

| method | run | unit | wall s | CPU-hours | peak RSS GiB | node |
|---|---|---|---|---|---|---|
| vsearch (n_mismatch) | SRR14351851 | R1 | 4424 | 18.341 | 142.9 | mjolnircomp01fl |
| vsearch (n_mismatch) | SRR14351851 | R2 | 4724 | 19.614 | 143.0 | mjolnircomp01fl |
| vsearch (n_mismatch) | ERR6001161 | R1 | 5178 | 20.643 | 141.2 | mjolnircomp06fl |
| vsearch (n_mismatch) | ERR6001161 | R2 | 5301 | 21.033 | 141.3 | mjolnircomp06fl |
| vsearch (n_mismatch) | ERR6001173 | R1 | 1616 | 5.927 | 141.2 | mjolnirgpu03fl |
| vsearch (n_mismatch) | ERR6001173 | R2 | 1668 | 6.161 | 141.2 | mjolnirgpu03fl |
| minimap2 preset | SRR14351851 | R1 | 102 | 0.072 | 32.6 | mjolnirgpu03fl |
| minimap2 preset | SRR14351851 | R2 | 98 | 0.069 | 31.6 | mjolnirgpu03fl |
| minimap2 preset | ERR6001161 | R1 | 68 | 0.049 | 33.3 | mjolnircomp15fl |
| minimap2 preset | ERR6001161 | R2 | 68 | 0.049 | 32.2 | mjolnircomp15fl |
| minimap2 preset | ERR6001173 | R1 | 101 | 0.071 | 31.8 | mjolnirgpu03fl |
| minimap2 preset | ERR6001173 | R2 | 100 | 0.071 | 31.4 | mjolnirgpu03fl |
| minimap2 sensitive | SRR14351851 | R1 | 83 | 0.059 | 32.7 | mjolnircomp12fl |
| minimap2 sensitive | SRR14351851 | R2 | 84 | 0.059 | 31.9 | mjolnircomp12fl |
| minimap2 sensitive | ERR6001161 | R1 | 83 | 0.057 | 32.4 | mjolnircomp01fl |
| minimap2 sensitive | ERR6001161 | R2 | 80 | 0.056 | 32.5 | mjolnircomp01fl |
| minimap2 sensitive | ERR6001173 | R1 | 135 | 0.097 | 32.5 | mjolnircomp06fl |
| minimap2 sensitive | ERR6001173 | R2 | 133 | 0.097 | 31.5 | mjolnircomp06fl |
| Kraken2 (conf 0.0) | SRR14351851 | pair | 2 | 0.001 | 1.3 | mjolnircomp15fl |
| Kraken2 (conf 0.0) | ERR6001161 | pair | 2 | 0.001 | 0.6 | mjolnircomp15fl |
| Kraken2 (conf 0.0) | ERR6001173 | pair | 3 | 0.001 | 0.6 | mjolnircomp01fl |
| Kraken2 (conf 0.1) | SRR14351851 | pair | 2 | 0.001 | 0.7 | mjolnircomp03fl |
| Kraken2 (conf 0.1) | ERR6001161 | pair | 1 | 0.001 | 0.6 | mjolnircomp03fl |
| Kraken2 (conf 0.1) | ERR6001173 | pair | 1 | 0.001 | 0.7 | mjolnircomp02fl |

BarcodeMapper and best-hit (sacct of the pilot jobs, **not comparable to classifier-only timings**: they include Snakemake, 8-CPU Bowtie2 mapping with a cap of 200 hits or `--all-hits`, and extraction; best-hit is a post-processing script on the BarcodeMapper BAM-derived extract, 1 CPU):

| pilot | mapping worker job | elapsed | MaxRSS (step) | alloc CPUs | controller elapsed |
|---|---|---|---|---|---|
| SRR14351851-pilot-standard | 45939594 | 00:02:50 | 24060477K | 8 | 00:11:04 |
| SRR14351851-pilot-allhits | 45939590 | 00:06:36 | 30937711K | 8 | 00:15:45 |
| ERR6001161-pilot-standard | 45939595 | 00:01:34 | 21973822K | 8 | 00:11:49 |
| ERR6001161-pilot-allhits | 45939581 | 02:01:35 | 33103791K | 8 | 02:10:17 |
| ERR6001173-pilot-standard | 45939591 | 00:01:31 | 21548698K | 8 | 00:11:09 |
| ERR6001173-pilot-allhits | 45939592 | 00:03:06 | 32663750K | 8 | 00:12:32 |

| best-hit run | job | state | elapsed | MaxRSS | alloc CPUs |
|---|---|---|---|---|---|
| SRR14351851 | 45939557 | COMPLETED | 00:00:02 | 796K | 1 |
| ERR6001161 | 45939562 | COMPLETED | 00:00:02 | 798K | 1 |
| ERR6001173 | 45939567 | COMPLETED | 00:00:02 | 841K | 1 |

## Consistency checks

- SRR14351851: best_hit hit set == BarcodeMapper standard: True; vsearch (n_mismatch) hit-only 220, BarcodeMapper-only 1; minimap2 preset hit-only 486, BarcodeMapper-only 1; Kraken2 0.0 hit-only 6869, BarcodeMapper-only 4; Kraken2 0.1 subset of 0.0: True
- ERR6001161: best_hit hit set == BarcodeMapper standard: True; vsearch (n_mismatch) hit-only 3, BarcodeMapper-only 0; minimap2 preset hit-only 57, BarcodeMapper-only 0; Kraken2 0.0 hit-only 703, BarcodeMapper-only 0; Kraken2 0.1 subset of 0.0: True
- ERR6001173: best_hit hit set == BarcodeMapper standard: True; vsearch (n_mismatch) hit-only 1, BarcodeMapper-only 1; minimap2 preset hit-only 16, BarcodeMapper-only 0; Kraken2 0.0 hit-only 1138, BarcodeMapper-only 2; Kraken2 0.1 subset of 0.0: True

## Exploratory notes (read-only, login node)

(a) Why does vsearch (n_mismatch) hit 446 fungal pairs against 227 for BarcodeMapper? SRR14351851: 220 pairs are hit by vsearch only and 1 by BarcodeMapper only. For 10 evenly spaced vsearch-only pairs (full listing: exploratory_notes.txt), the vsearch alignments are partial: 81 to 214 aligned bases, 27 to 71% of the 300 bp read, at 97 to 100% identity. BarcodeMapper maps with Bowtie2 in end-to-end mode (no --local), so the whole read has to align. In 4 of the 10 pairs Bowtie2 reported no alignment at all; in the other 6 it aligned the full read (about 300 bases) to other references with NM 11 to 35 (3.7 to 11.7% of covered bases), which fails the 3% rule. So the extra vsearch hits come from local, partial alignments of 300 bp reads that overhang the marker; the shared absolute 80-base threshold admits them, and this is not a -k 200, soft-clipping or threshold effect in Bowtie2. The same logic probably applies to minimap2 (local alignments); this was not checked for minimap2 or the animal mock.

(b) Why does Kraken2 classify few pairs? Pairs with at least one 35-mer present in the database (a non-zero, non-ambiguous entry in the k-mer hit string of output.txt), against pairs classified: SRR14351851 13,818 (15.0%) vs 7,151 (7.7%); ERR6001161 9,146 (9.3%) vs 757 (0.8%); ERR6001173 8,247 (8.3%) vs 1,162 (1.2%). Only 0.74%, 0.17% and 0.13% of read k-mers hit the database. So most pairs have no marker k-mer (shotgun reads are mostly off-target, as for BarcodeMapper), and many pairs with a k-mer hit stay unclassified; the command does not override Kraken2's default --minimum-hit-groups 2, which is the probable cause (not verified). At confidence 0.1 classified pairs fall to 1,096, 96 and 26. Kraken2 is not low relative to BarcodeMapper in the animal mock and blank (757 and 1,162 classified vs 21 and 4 pairs hit); those extra calls are counted in the tables above.

## Caveats

- One 100k-pair subsample per mock (92k to 99k pairs after fastp); counts of marker pairs are tiny (a few to a few hundred), so detection numbers are indicative only; no replicates or confidence intervals.
- Expected-taxa matching uses exact paper binomials; species-name synonyms (for example Nakaseomyces glabratus for Candida glabrata) count as unexpected, so recall is a lower bound.
- Matching ignores kingdom: off-target-kingdom calls count as unexpected. The `all-fungi` reference also contains Animalia records.
- vsearch is not exhaustive (`--maxrejects 1000`, `--maxaccepts 0`, top-identity ties only) and its alignments are global with free terminal gaps; the acceptance rule (covered >= 80 bases, NM <= 3%) is applied post hoc.
- minimap2 returns many more hit pairs than BarcodeMapper on the animal mock (78 vs 21) without a mechanism being established here.
- Kraken2 classifies few pairs; it needs at least 2 hit groups by default (see note b), and the database hash is only about 0.4 GB.
- BLASTn is pending. Synthetic per-fragment truth (Project 02) is not built; synthetic mode was not run.
- Blank ERR6001173 uses a placeholder expected species so that every call is unexpected.

## Appendix: invalid vsearch run (N-run artifact), not part of the comparison

The first vsearch run scored N as match; 89 to 95% of its accepted alignments lie on reference stretches that are at least 50% N (see evaluation/v2/n_run_diagnostic.tsv). Retained for the record only.

SRR14351851:

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| vsearch (invalid, N artifact) | 36889 | 36795 | 0.39829 | 12/38 (9) | 14/24 (12) | 12/14 (11) | 826 (27866) | 659 | 292 |

ERR6001161:

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| vsearch (invalid, N artifact) | 80786 | 80785 | 0.82303 | 5/28 (2) | 8/24 (5) | 7/18 (4) | 1277 (27818) | 983 | 371 |

ERR6001173:

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| vsearch (invalid, N artifact) | 48588 | 48588 | 0.49138 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 1505 (21906) | 1196 | 448 |

