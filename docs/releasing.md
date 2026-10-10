# Release procedure

1. Update `project.version` in `pyproject.toml` using `MAJOR.MINOR.PATCH`. Set the same version in `CITATION.cff`. Move any changes from `Unreleased` into that version's changelog section and date the section, leaving an empty `Unreleased` section.
2. Run the tests in the pinned environment, including the end-to-end checks: `python -m pytest -q`. Review the limitations and database provenance notes for the release.
3. Check the proposed tag and build artifacts:

   ```sh
   python scripts/check_release.py --tag vX.Y.Z
   python -m pip install build twine
   python -m build
   python -m twine check dist/*
   ```

4. Inspect the source archive for `LICENSE` and `CITATION.cff`, and the wheel for `LICENSE`, workflow files, profiles, and scripts. Install the wheel in a clean environment and check `barcodemapper --version` and `--help`.
5. After review, create the matching Git tag `vX.Y.Z`, publish a GitHub release from that tag, and upload the built wheel and source distribution. Do not bundle BOLD/UNITE derived data without confirming redistribution rights.
6. To publish on Bioconda, update `packaging/bioconda/meta.yaml` to the release version and SHA-256 of the uploaded source distribution, then submit it as `recipes/barcodemapper/meta.yaml` in [bioconda-recipes](https://github.com/bioconda/bioconda-recipes). Verify the channel package after the recipe is merged before advertising `conda install -c bioconda barcodemapper`. The recipe includes Snakemake, fastp, Bowtie2, and samtools as runtime dependencies.

`scripts/check_release.py` rejects a proposed tag if its version, citation version, or dated changelog section is inconsistent.
