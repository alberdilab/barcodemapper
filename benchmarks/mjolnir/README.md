# BarcodeMapper publication benchmarks on Mjolnir

These are four distinct benchmark projects. The projects share an immutable
software release and reference snapshot, while each has its own Conda prefix,
input manifest, working directory, results, and run record.

[Execution record](EXECUTION.md) tracks the Mjolnir jobs and frozen input metadata.

| Project | Plan | Main result |
| --- | --- | --- |
| 01 | [Defined shotgun mock communities](01-mock-communities.md) | Taxon detection and false calls on biological mixtures |
| 02 | [Reference novelty and ambiguity](02-reference-novelty.md) | Behavior when source taxa are absent or barcodes conflict |
| 03 | [Method comparison](03-method-comparison.md) | Accuracy and compute tradeoffs against controlled baselines |
| 04 | [Scaling](04-scaling.md) | Time, CPU, memory, and scratch use as data grow |

## Shared execution contract

1. Use `/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks` as
   `BENCH_ROOT`, or replace that path consistently in every command and record.
   Do not use the shared lab Conda environments for these runs.
2. Project 03 uses BarcodeMapper **v1.2.3** (see
   [03-method-comparison.md](03-method-comparison.md)); the v1.1.1 text below
   is the frozen history of the earlier projects. Freeze the BarcodeMapper **v1.1.1** wheel before starting. It was built
   from the working tree on top of Git commit
   `893040ff24eb85a706ce2b4bf97d12e55619a11d`; v1.1.0 (commit
   `7efbdd4941a74ac3b7b5e24ca420be6ba075aed1`) was superseded because its
   Slurm profile capped `index_db` at 2 threads and under-sized its memory.
   The local release wheel `barcodemapper-1.1.1-py3-none-any.whl` has SHA-256
   `67f31651097153ba7310ae26e08c3486876635e303d5c19b7b7386d56c967a67`.
   Verify the transferred or downloaded wheel with [wheel.sha256](wheel.sha256).
   If the published wheel differs, record the new digest and investigate before using
   it. Install the checked wheel with `python -m pip install --no-deps`.
3. Copy [environment.yml](environment.yml) into each project's `software/`
   directory. Create **one Conda prefix per project** from that file, then save
   `conda list --explicit` as `software/conda-linux-64.explicit.txt` and
   `conda list --json` as `software/conda-packages.json`. The explicit file is
   the replay lock for Linux x86_64. Preserve the environment YAML, the lock,
   the wheel, and the wheel checksum with the results.
4. Download BOLD Public and UNITE All Eukaryotes source FASTAs **once**, save
   them under `BENCH_ROOT/references/sources/`, and record release identifiers,
   DOIs, source URLs, and SHA-256 hashes. `barcodemapper --build` without `-b`
   and `-u` chooses whatever is current on the day of execution; for the
   benchmark, build from the retained source snapshot with explicit `-b` and
   `-u` paths. The current BOLD release has zero-length records that fail the
   v1.1.0 reference validator. Run
   [drop_empty_source_records.py](drop_empty_source_records.py) once and use its
   read-only, checksum-traced FASTAs in `references/cleaned-sources/`. The
   original downloads remain unchanged; the filter removes only records with
   no sequence bases. Keep the resulting FASTA, all six Bowtie2 index files, taxonomy
   overview, and build manifest read-only. The repository's
   [freeze_sources.py](freeze_sources.py) acquires the source FASTAs and writes
   `source-releases.json` and `SHA256SUMS`. Run discovery only once; later runs
   use the retained FASTAs, never another "latest" lookup.
5. Build two clearly named database snapshots from those same sources:
   `default.fa.gz` with the package's default BOLD/UNITE taxon filters, and
   `all-fungi.fa.gz` with `-x k__Animalia -y k__Fungi`. The latter is needed for
   the broad fungal mock community; many of its yeasts fall outside the
   package's default ITS classes. Record both FASTA SHA-256 hashes. Do not
   compare methods run against different snapshots.
6. Save one filled [run record](run-record.template.yaml) per invocation. Use
   absolute input paths, the input checksums, database checksum, wheel hash,
   Conda lock hash, random seed, full command, Slurm job IDs, and output hashes.
   Preserve BarcodeMapper's own `.manifest.json` beside each result.

### Environment setup on Mjolnir

The `useconda` alias is only available in an interactive shell. From an
interactive Mjolnir session, the setup for each project is shown below.
Submit this work through [setup_environment.sbatch](setup_environment.sbatch)
for Project 01; [clone_environment.sbatch](clone_environment.sbatch) replays
the resulting explicit lock for Projects 02–04 on compute nodes:

```bash
BENCH_ROOT=/projects/alberdilab/scratch/jpl786/projects/barcodemapper-benchmarks
PROJECT=01-mock-communities  # replace for each project
PROJECT_ROOT="$BENCH_ROOT/projects/$PROJECT"
mkdir -p "$PROJECT_ROOT"/{software,inputs,work,results,logs}
cp "$BENCH_ROOT/plans/environment.yml" "$PROJECT_ROOT/software/environment.yml"
(cd "$BENCH_ROOT/software" && sha256sum -c "$BENCH_ROOT/plans/wheel.sha256")
useconda
conda env create --prefix "$PROJECT_ROOT/software/env" --file "$PROJECT_ROOT/software/environment.yml"
conda activate "$PROJECT_ROOT/software/env"
python -m pip install --no-deps "$BENCH_ROOT/software/barcodemapper-1.1.1-py3-none-any.whl"
conda list --explicit > "$PROJECT_ROOT/software/conda-linux-64.explicit.txt"
conda list --json > "$PROJECT_ROOT/software/conda-packages.json"
barcodemapper --version
```

To replay an environment after the first solve, create a new prefix with
`conda create --prefix NEW_PREFIX --file conda-linux-64.explicit.txt`, then
install the same checked wheel with `python -m pip install --no-deps`. Keep
the explicit lock and wheel together; the lock alone does not include the pip
installation.

Transfer this entire `benchmarks/mjolnir/` directory to `BENCH_ROOT/plans/`
and the verified wheel to `BENCH_ROOT/software/` before using the commands.
For a non-interactive SSH or Slurm shell, activate with
`bash -ic 'useconda; conda activate /absolute/prefix; ...'`. Do not place
credentials in the project or logs.

Mjolnir's compute-node Slurm configuration currently prevents `sacct` queries
from the pinned Snakemake executor. The tested [accounting wrapper](sacct-wrapper.sh)
and its captured site configuration are described in the
[execution record](EXECUTION.md). Link the wrapper into each benchmark
environment's `bin/sacct` before submitting Snakemake controllers.

For long analyses, submit the bundled [controller.sbatch](controller.sbatch)
instead of leaving a Snakemake controller attached to a login shell. Export
the project's `BENCH_ENV` first, then pass the exact BarcodeMapper command as
arguments. Site-specific account and partition flags belong in the `sbatch`
invocation and must be saved in the run record:

```bash
export BENCH_ENV="$PROJECT_ROOT/software/env"
sbatch --export=ALL --job-name=bm-mock-animal01 \
  "$BENCH_ROOT/plans/controller.sbatch" \
  barcodemapper --slurm \
  -1 "$PROJECT_ROOT/inputs/pairs/animal01_1.fastq.gz" \
  -2 "$PROJECT_ROOT/inputs/pairs/animal01_2.fastq.gz" \
  -d "$BENCH_ROOT/references/default.fa.gz" \
  -o "$PROJECT_ROOT/results/animal01-standard.tsv" \
  -t "$PROJECT_ROOT/work/animal01-standard" \
  --standard --adapter-r1 auto --adapter-r2 auto
```

### Database build, performed once

Acquire and verify the source FASTAs with
`sbatch --export=ALL,BENCH_ROOT="$BENCH_ROOT",BENCH_ENV="$BENCH_ENV" "$BENCH_ROOT/plans/freeze_sources.sbatch"`.
The job downloads and decompresses the sources on `cpuqueue`. Once it succeeds,
run `drop_empty_source_records.sbatch`; verify its
`empty-record-filter.json` and exclusion count. Then build with fixed filters.
Submit each long `barcodemapper --build --slurm`
invocation through `controller.sbatch` as above, passing the build commands
below as its arguments:

```bash
cd "$BENCH_ROOT/references/sources"
sha256sum -c SHA256SUMS

barcodemapper --build --slurm \
  -b "$BENCH_ROOT/references/cleaned-sources/bold.fasta" \
  -u "$BENCH_ROOT/references/cleaned-sources/unite.fasta" \
  -d "$BENCH_ROOT/references/default.fa.gz" \
  -t "$BENCH_ROOT/references/work-default"

barcodemapper --build --slurm \
  -b "$BENCH_ROOT/references/cleaned-sources/bold.fasta" \
  -u "$BENCH_ROOT/references/cleaned-sources/unite.fasta" \
  -x k__Animalia -y k__Fungi \
  -d "$BENCH_ROOT/references/all-fungi.fa.gz" \
  -t "$BENCH_ROOT/references/work-all-fungi"
```

Run a small paired-end pilot before the full projects. The packaged Slurm
profile caps rule threads at eight and may need site-specific resource
adjustment before large runs; record any such adjustment as a versioned
project artifact. Preserve the original profile and the adjusted copy.

## Interpretation rules

- Count each read **pair/fragment** once. fastp's `total_reads` counts mates,
  so divide it by two only after confirming both mates survived as a pair.
- BarcodeMapper's table propagates counts to ancestors. Never sum all rows to
  infer a total. Use per-fragment assignment files under `work/taxonomy/` for
  read-level evaluation.
- A correct genus assignment for a species of known truth is a genus success
  and a species abstention, not a correct species call. Report resolution and
  wrong overcalls separately from rank-specific precision and recall.
- Biological mock-community membership is known; the exact origin of every
  shotgun read is generally not. Use the in-silico project for per-read truth.
- Do not describe marker-matching fragment counts as organism abundance or
  biomass. Any count-versus-input-mass association is secondary.
