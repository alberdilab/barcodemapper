# Testing BarcodeMapper

From the repository root, in an environment with the package and `pytest` installed:

```sh
python -m pytest -q
```

The unit tests check input validation, per-alignment edit filtering, taxonomy aggregation, zero-hit samples, report accounting, and handling of conflicting BOLD records. Two end-to-end tests generate a tiny deterministic paired-read fixture and use real Snakemake, fastp, Bowtie2, and samtools executables. They check that an ambiguous species match stops at genus, an unrelated sample remains at zero, a combined database can be built, and source checksums are recorded. Integration tests skip when required executables are absent. Run only unit tests with `python -m pytest -q -m 'not integration'`.

These checks verify software behavior on synthetic data. They do not establish sensitivity, specificity, or abundance accuracy on biological samples; see [validation scope](docs/validation.md).
