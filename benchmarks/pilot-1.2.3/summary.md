# Pilots 1.2.3: comparator evaluation on real mocks (exploratory)

BarcodeMapper 1.2.3 `--standard` (-k 200), `--all-hits`, single best-hit ablation, vsearch, minimap2 preset (primary minimap2 row) and minimap2 sensitive variant (`-f 20000,100000`, separately labelled). Truth: Project 01 `expected-taxa.tsv` (exact paper binomials, synonyms unreviewed, so recall is a lower bound and synonym names count as unexpected). Matching is by name at species, genus and family, regardless of kingdom (so off-target kingdoms count as unexpected). Fragments = read pairs. Counts are one pilot per marker, 92k to 99k pairs, so they are indicative only. Per-fragment synthetic truth (Project 02) is not built; synthetic mode was not run. **Pending: Kraken2 (database builds queued) and BLASTn (subset awaiting decisions).**

## SRR14351851: fungal pooled DNA, all-fungi reference

Detected / expected (min 1 fragment; in parentheses min 2 fragments); unexpected = taxa called that are not expected, with total fragments.

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 227 | 181 | 0.00196 | 5/38 (2) | 9/24 (4) | 10/14 (4) | 8 (88) | 8 | 6 |
| barcodemapper_allhits | 226 | 180 | 0.00195 | 5/38 (2) | 9/24 (3) | 10/14 (4) | 7 (87) | 8 | 6 |
| best_hit | 227 | 227 | 0.00246 | 11/38 (7) | 13/24 (9) | 12/14 (8) | 17 (115) | 9 | 7 |
| vsearch | 36889 | 36795 | 0.39829 | 12/38 (9) | 14/24 (12) | 12/14 (11) | 826 (27866) | 659 | 292 |
| minimap2_preset | 712 | 577 | 0.00625 | 3/38 (2) | 4/24 (2) | 6/14 (4) | 18 (186) | 25 | 23 |
| minimap2_sensitive | 716 | 569 | 0.00616 | 3/38 (2) | 4/24 (2) | 6/14 (4) | 18 (185) | 25 | 23 |

Top unexpected species per method (fragments):

- barcodemapper_standard: Zygotorulaspora_dagestanica (58); Cordyceps_fumosorosea (17); Nakaseomyces_glabratus (3); Scedosporium_prolificans (3); Vishniacozyma_pseudopenaeus (3)
- barcodemapper_allhits: Zygotorulaspora_dagestanica (58); Cordyceps_fumosorosea (17); Nakaseomyces_glabratus (3); Scedosporium_prolificans (3); Vishniacozyma_pseudopenaeus (3)
- best_hit: Zygotorulaspora_dagestanica (59); Cordyceps_fumosorosea (18); Vishniacozyma_pseudopenaeus (8); Holtermannia_saccardoi (4); Nakaseomyces_glabratus (4)
- vsearch: Tubastraea_coccinea (10680); Ceratosolen_jucundus (1107); Ravinia_derelicta (434); Rana_temporaria (427); Muhabbetana_insulata (422)
- minimap2_preset: Zygotorulaspora_dagestanica (125); Cordyceps_fumosorosea (21); Chrysomycena_dunicola (7); Holtermannia_saccardoi (7); Scutellospora_castanea (6)
- minimap2_sensitive: Zygotorulaspora_dagestanica (125); Cordyceps_fumosorosea (20); Chrysomycena_dunicola (7); Holtermannia_saccardoi (7); Scutellospora_castanea (6)

## ERR6001161: animal mock, default reference

Detected / expected (min 1 fragment; in parentheses min 2 fragments); unexpected = taxa called that are not expected, with total fragments.

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 21 | 21 | 0.00021 | 5/28 (2) | 8/24 (4) | 8/18 (3) | 3 (6) | 2 | 4 |
| barcodemapper_allhits | 21 | 21 | 0.00021 | 5/28 (2) | 8/24 (4) | 8/18 (3) | 3 (6) | 2 | 4 |
| best_hit | 21 | 21 | 0.00021 | 5/28 (2) | 9/24 (4) | 9/18 (3) | 6 (10) | 2 | 4 |
| vsearch | 80786 | 80785 | 0.82303 | 5/28 (2) | 8/24 (5) | 7/18 (4) | 1277 (27818) | 983 | 371 |
| minimap2_preset | 78 | 74 | 0.00075 | 4/28 (2) | 7/24 (4) | 8/18 (4) | 4 (20) | 4 | 7 |
| minimap2_sensitive | 79 | 75 | 0.00076 | 4/28 (2) | 7/24 (4) | 8/18 (4) | 4 (20) | 4 | 8 |

Top unexpected species per method (fragments):

- barcodemapper_standard: Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Steganacarus_spinosus (1)
- barcodemapper_allhits: Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Steganacarus_spinosus (1)
- best_hit: Eupelops_sp._SFH2024_01 (3); Poikilolaimus_oxycercus (2); Seinura_italiensis (2); Enchytraeus_bigeminus_L1 (1); Paramacrobiotus_fairbanksi (1)
- vsearch: Ceratosolen_jucundus (2292); Tubastraea_coccinea (2188); Chrysis_solida (998); Idarnes_dimorphicus (573); Cecidonius_sp. (415)
- minimap2_preset: Hemiberlesia_rapax (14); Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Enchytraeus_sp._CEF (1)
- minimap2_sensitive: Hemiberlesia_rapax (14); Eupelops_sp._SFH2024_01 (3); Seinura_italiensis (2); Enchytraeus_sp._CEF (1)

## ERR6001173: extraction blank (negative control), default reference

Detected / expected (min 1 fragment; in parentheses min 2 fragments); unexpected = taxa called that are not expected, with total fragments.

| method | pairs with hit | classified | classified fraction | species detected | genus detected | family detected | unexpected species (fragments) | unexpected genera | unexpected families |
|---|---|---|---|---|---|---|---|---|---|
| barcodemapper_standard | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |
| barcodemapper_allhits | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |
| best_hit | 4 | 4 | 0.00004 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 2 (2) | 2 | 2 |
| vsearch | 48588 | 48588 | 0.49138 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 1505 (21906) | 1196 | 448 |
| minimap2_preset | 20 | 18 | 0.00018 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 1 (1) | 1 | 2 |
| minimap2_sensitive | 21 | 19 | 0.00019 | 0/1 (0) | 0/1 (0) | 0/0 (0) | 0 (0) | 0 | 1 |

Top unexpected species per method (fragments):

- barcodemapper_standard: none
- barcodemapper_allhits: none
- best_hit: Deuterostichococcus_epilithicus (1); Sphaeridia_pumilis (1)
- vsearch: Tubastraea_coccinea (1805); Ceratosolen_jucundus (830); Rana_temporaria (426); Chrysis_solida (274); Homo_sapiens (231)
- minimap2_preset: Chironomus_riparius (1)
- minimap2_sensitive: none

## Exploratory timings (not the final timing runs)

16 CPUs, classifier-only (index/database load included), `time_classifier.sh`; nodes vary (including a GPU node), so wall times are not strictly comparable. BarcodeMapper, best-hit and Kraken2 timings are not part of this table.

| method | run | mate | wall s | CPU-hours | peak RSS GiB | node |
|---|---|---|---|---|---|---|
| vsearch | SRR14351851 | R1 | 5848 | 24.31 | 142.9 | mjolnirgpu03fl |
| vsearch | SRR14351851 | R2 | 5624 | 23.49 | 142.9 | mjolnirgpu03fl |
| vsearch | ERR6001161 | R1 | 2342 | 9.73 | 141.2 | mjolnircomp14fl |
| vsearch | ERR6001161 | R2 | 2351 | 9.77 | 141.3 | mjolnircomp14fl |
| vsearch | ERR6001173 | R1 | 1199 | 4.60 | 141.3 | mjolnircomp15fl |
| vsearch | ERR6001173 | R2 | 1196 | 4.60 | 141.2 | mjolnircomp15fl |
| minimap2 preset | SRR14351851 | R1 | 102 | 0.07 | 32.6 | mjolnirgpu03fl |
| minimap2 preset | SRR14351851 | R2 | 98 | 0.07 | 31.6 | mjolnirgpu03fl |
| minimap2 preset | ERR6001161 | R1 | 68 | 0.05 | 33.3 | mjolnircomp15fl |
| minimap2 preset | ERR6001161 | R2 | 68 | 0.05 | 32.2 | mjolnircomp15fl |
| minimap2 preset | ERR6001173 | R1 | 101 | 0.07 | 31.8 | mjolnirgpu03fl |
| minimap2 preset | ERR6001173 | R2 | 100 | 0.07 | 31.4 | mjolnirgpu03fl |
| minimap2 sensitive | SRR14351851 | R1 | 83 | 0.06 | 32.7 | mjolnircomp12fl |
| minimap2 sensitive | SRR14351851 | R2 | 84 | 0.06 | 31.9 | mjolnircomp12fl |
| minimap2 sensitive | ERR6001161 | R1 | 83 | 0.06 | 32.4 | mjolnircomp01fl |
| minimap2 sensitive | ERR6001161 | R2 | 80 | 0.06 | 32.5 | mjolnircomp01fl |
| minimap2 sensitive | ERR6001173 | R1 | 135 | 0.10 | 32.5 | mjolnircomp06fl |
| minimap2 sensitive | ERR6001173 | R2 | 133 | 0.10 | 31.5 | mjolnircomp06fl |

## Consistency checks (for converter or evaluator bugs)

- SRR14351851: best_hit hit pairs = 227, barcodemapper_standard = 227, identical hit sets: True; allhits-only pairs 0, standard-only pairs 1
- SRR14351851: pairs hit by vsearch and not by BarcodeMapper standard: 36662; BarcodeMapper standard hit and not vsearch: 0; minimap2 preset hit and not standard: 486; minimap2 sensitive hit and not standard: 490
- SRR14351851: minimap2 preset hit pairs 712, sensitive 716, preset not in sensitive 0
- ERR6001161: best_hit hit pairs = 21, barcodemapper_standard = 21, identical hit sets: True; allhits-only pairs 0, standard-only pairs 0
- ERR6001161: pairs hit by vsearch and not by BarcodeMapper standard: 80765; BarcodeMapper standard hit and not vsearch: 0; minimap2 preset hit and not standard: 57; minimap2 sensitive hit and not standard: 58
- ERR6001161: minimap2 preset hit pairs 78, sensitive 79, preset not in sensitive 0
- ERR6001173: best_hit hit pairs = 4, barcodemapper_standard = 4, identical hit sets: True; allhits-only pairs 0, standard-only pairs 0
- ERR6001173: pairs hit by vsearch and not by BarcodeMapper standard: 48585; BarcodeMapper standard hit and not vsearch: 1; minimap2 preset hit and not standard: 16; minimap2 sensitive hit and not standard: 17
- ERR6001173: minimap2 preset hit pairs 20, sensitive 21, preset not in sensitive 0

## Notes

- Blank ERR6001173 has no expected taxa; evaluate_calls.py needs at least one row, so a placeholder species `NO_EXPECTED_TAXA_BLANK_CONTROL` was used. All calls are therefore unexpected.
- vsearch uses `--maxrejects 1000` (not exhaustive) and the exact BarcodeMapper acceptance rule applied post hoc; minimap2 rows use the same rule.
- All conversions here use the committed (HEAD) converter; `script_and_truth_hashes.txt` lists hashes.

## WARNING: the vsearch row is dominated by an N-run artifact (not a real difference)

vsearch reports 40 to 82% of pairs as hits (BarcodeMapper standard: 0.02 to 0.25%), with thousands of unexpected species in all three samples including the blank (Homo sapiens, Rana temporaria, Tubastraea coccinea). Diagnostic on ERR6001161 R1: 111,992 of 118,461 accepted hits (94.5%) align onto reference stretches that contain N (all of them at least 50% N). vsearch scores N as a match, so a read aligns at 100% identity to an N run, and the exact NM filter applied afterwards cannot see it because the N columns are not in vsearch's mismatch count. With --top_hits_only these 100%-identity N hits also hide the real 97 to 99% hits. The vsearch row should be treated as invalid. Fix (needs new jobs, not run): rerun vsearch with --n_mismatch (N counted as mismatch) or an N-masked/N-free reference copy, as a new frozen config. minimap2 and BarcodeMapper rows are not affected by this mechanism. Detail: vsearch_N_artifact_diagnostic.txt.
