# Release procedure

1. Update `project.version` in `pyproject.toml` using `MAJOR.MINOR.PATCH`. Set the same version in `CITATION.cff`. Move any changes from `Unreleased` into that version's changelog section and date the section, leaving an empty `Unreleased` section.
2. Run the tests in the pinned environment, including the end-to-end checks: `python -m pytest -q`. Review the limitations and database provenance notes for the release.
3. Check whether newer dependencies exist and whether results survive them. From a machine with the Mjolnir connection:

   ```sh
   scripts/release_dependency_check.sh
   ```

   It compares the pins in `barcodemapper_env.yml` with the latest conda-forge/bioconda releases. If anything is newer, it builds a throwaway environment with the latest versions on Mjolnir, installs the checkout, runs the small regression dataset (`tests/regression/`) and requires the taxonomy table to be identical to the committed golden table. It repeats the check through Slurm (`--no-slurm` skips that) and `--record` appends the outcome to [](dependency-checks.md), which also lists every check and its command. A pass means the newer versions can be adopted by editing the pins; a failure means they must wait until the difference is understood. Run it with `--record` and commit the new row. `--force` runs the HPC check even when every pin is current, and `BM_HPC_SSH` sets the ssh command if `ssh mjolnirgate.unicph.domain` is not enough. The same check runs weekly on GitHub Actions (`dependency-check.yml`).
4. Check the proposed tag and build artifacts:

   ```sh
   python scripts/check_release.py --tag vX.Y.Z
   python -m pip install build twine
   python -m build
   python -m twine check dist/*
   ```

5. Inspect the source archive for `LICENSE` and `CITATION.cff`, and the wheel for `LICENSE`, workflow files, profiles, and scripts. Install the wheel in a clean environment and check `barcodemapper --version` and `--help`.
6. After review, create the matching Git tag `vX.Y.Z`, publish a GitHub release from that tag, and upload the built wheel and source distribution. Do not attach BOLD/UNITE derived data to the GitHub release; databases are published on Zenodo (see below).
7. To publish on Bioconda, update `packaging/bioconda/meta.yaml` to the release version and SHA-256 of the uploaded source distribution, then submit it as `recipes/barcodemapper/meta.yaml` in [bioconda-recipes](https://github.com/bioconda/bioconda-recipes). Verify the channel package after the recipe is merged before advertising `conda install -c bioconda barcodemapper`. The recipe includes Snakemake, fastp, Bowtie2, and samtools as runtime dependencies.

`scripts/check_release.py` rejects a proposed tag if its version, citation version, or dated changelog section is inconsistent.

## Publishing a database release to Zenodo

Database releases are independent of software releases. See `scripts/README.md` for the full options.

1. On Mjolnir, build and upload a draft to the Zenodo sandbox first, using `ZENODO_TOKEN` from `~/.tokens` (never on a command line):

   ```sh
   sbatch --export=NONE scripts/build_and_publish.sh --db-dir /projects/alberdilab/scratch/jpl786/barcodemapper-dbs \
     --env /projects/alberdilab/data/environments/conda/barcodemapper --sandbox
   ```

2. Inspect the draft in Zenodo's UI: files, checksums, licence (CC BY-SA 4.0), and the attribution to BOLD and UNITE.
3. Repeat without `--sandbox` and publish with `--publish`, or publish the production draft in the UI. Zenodo drafts block new versions until they are published or deleted, so clean up unused drafts.
4. Set `DEFAULT_RECORD` in `barcodemapper/zenodo.py` to the printed concept record id, so that `--fetch` finds the latest version by default.
5. Add the record to the README and `CITATION.cff` as appropriate, then follow the release procedure above.
