# Development

## Tests

From the repository root, in an environment with the package and `pytest` installed:

```sh
python -m pytest -q
```

The unit tests check input validation, reference header validation, exact read ID handling, per-alignment edit filtering, taxonomy aggregation, zero-hit samples, report accounting, and handling of conflicting BOLD records. The integration checks exercise two end-to-end workflows with real Snakemake, fastp, Bowtie2, and samtools executables. They check that an ambiguous species match stops at genus, an unrelated sample remains at zero, a combined database can be built, and source checksums are recorded. Integration tests skip when required executables are absent. Run only unit tests with `python -m pytest -q -m 'not integration'`.

These checks verify software behavior on synthetic data. They do not establish sensitivity, specificity, or abundance accuracy on biological samples; see [](limitations.md).

## Dependency and regression checks

`tests/regression/` holds a small synthetic dataset and the taxonomy table it must produce. `scripts/regression_check.py` runs the installed pipeline on it and fails on any difference, and prints the tool versions used:

```sh
python scripts/regression_check.py
```

`scripts/check_dependencies.py` lists the pins in `barcodemapper_env.yml` next to the latest released versions and can write an environment file with the latest versions (`--write-candidate`). `scripts/release_dependency_check.sh` combines the two on Mjolnir; see [](releasing.md). A scheduled GitHub Actions workflow repeats the check weekly. Adopting new pins stays a manual edit of `barcodemapper_env.yml` after a passing check.

## Building these docs

```sh
python -m pip install -r docs/requirements.txt
python -m sphinx -W -b html docs docs/_build/html
```

## Releases

```{include} releasing.md
:start-line: 2
```
