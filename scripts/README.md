# Publishing a database to Zenodo

Maintainers run `python scripts/publish_zenodo.py` after `barcodemapper --build`.
The script stages the files expected by `barcodemapper --fetch`, verifies their
checksums, and creates a Zenodo draft. It publishes only with `--publish`.

The derived database is distributed under CC BY-SA 4.0, with attribution to
the BOLD (CC BY-SA 4.0) and UNITE (CC BY 4.0) sources in `ATTRIBUTION.txt`.
Review source licences and release IDs before publishing. Automatic source
provenance is required.

When calling `publish_zenodo.py` directly, set `ZENODO_TOKEN` in its environment.
For a first upload, use a sandbox token and inspect the plan and the resulting draft:

```sh
python scripts/publish_zenodo.py --database /path/bmdb_20261021 --sandbox --dry-run
python scripts/publish_zenodo.py --database /path/bmdb_20261021 --sandbox
```

Inspect and publish a draft in Zenodo's UI. Alternatively, `--publish` creates
and publishes a new draft in one run. Pass `--record CONCEPT_ID` to create a new
version of an existing deposition. The published output prints the concept
record ID to put into `barcodemapper/zenodo.py` as `DEFAULT_RECORD` so `--fetch`
can find the latest version by default. Use `--include-index` to upload the six
Bowtie2 index files.

## Build and publish on Mjolnir

`scripts/build_and_publish.sh` builds the current BOLD+UNITE database in scratch
and uploads it as a Zenodo draft. Use a pinned BarcodeMapper Conda environment
with Snakemake and Bowtie2; `--env` accepts its name or absolute prefix under
`/projects/alberdilab/data/environments/conda/`. Set `ZENODO_SANDBOX` for
`--sandbox` and `ZENODO_PRODUCTION` otherwise in Codeome Settings > Secrets,
enabled for HPC and this project. Mjolnir jobs use `sbatch --export=NONE` and
load secrets inside the job with `source ~/.codeome/secrets/<project>.env`;
this project's file is `~/.codeome/secrets/barcodemapper-91de8232.env`.
The wrapper reads only the selected variable from the environment or that file,
falling back to `~/.tokens`. Set `BM_SECRETS_FILE` to try a different secrets
file first. It does not accept a bare `ZENODO_TOKEN`. Never put a token on a
command line. Check that compute nodes can reach BOLD, UNITE, and Zenodo.

```sh
sbatch --export=NONE scripts/build_and_publish.sh --db-dir /projects/alberdilab/scratch/jpl786/barcodemapper-dbs \
  --env /projects/alberdilab/data/environments/conda/barcodemapper --sandbox
```

The sbatch header starts at 8 CPUs, 256 GB RAM, and three days. Adjust it for
the source size and queue limits: Bowtie2 indexing needs roughly 28 times the
compressed FASTA size in RAM; allow disk for source FASTAs, the combined FASTA,
index, and upload staging. `--reuse /path/to/bmdb_<date>` skips the build.
`--dry-run` requires an existing database, runs the publisher's offline checks,
and does not need a token. Review the draft before publishing it in Zenodo's UI;
`--publish` instead creates and immediately publishes a new draft. Afterward,
record the printed concept ID in `barcodemapper/zenodo.py` as `DEFAULT_RECORD`.
