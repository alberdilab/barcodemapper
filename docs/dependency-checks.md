# Dependency checks

Every release checks whether newer versions of the pinned dependencies exist and whether the pipeline still gives identical results with them. The checks and how to repeat them are described under Releases in [](development.md); this page records each run.

## Checks

| Check | Where | Command | Passes when |
| --- | --- | --- | --- |
| Unit and integration tests | local or CI | `python -m pytest -q` | all tests pass (integration tests skip without the workflow executables) |
| Version lookup | local | `python scripts/check_dependencies.py` | lists pinned versus latest; exit 0 |
| Regression, local profile | Mjolnir | `scripts/release_dependency_check.sh` | taxonomy table identical to `tests/regression/expected_taxonomy.tsv` |
| Regression, Slurm profile | Mjolnir | same command (add `--no-slurm` to skip) | same, with the jobs submitted through Slurm |

`scripts/release_dependency_check.sh --record` appends a row below. Run it with `--force` if every pin is already current, and `BM_HPC_SSH` if `ssh mjolnirgate.unicph.domain` is not enough. The golden table changes only on purpose: regenerate it with `scripts/regression_check.py --update-golden` in the pinned environment and explain why in the changelog.

## Results

| Date | BarcodeMapper (commit) | Local profile | Slurm profile | Versions tested | Outcome |
| --- | --- | --- | --- | --- | --- |
| 2026-10-18 | 1.2.2 (7a114c3) | PASS | PASS | biopython 1.88, bowtie2 2.5.5, fastp 1.4.0, numpy 2.5.4, pandas 2.3.3, plotly 7.1.0, pysam 0.24.1, python 3.12.15, pyyaml 6.0.3, samtools 1.24, snakemake 9.27.0, snakemake-executor-plugin-slurm 2.8.0 | PASS |
