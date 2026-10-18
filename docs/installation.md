# Installation

Python 3.12 or newer is required.

## Released wheel

Create an environment with the workflow tools and install the wheel:

```sh
conda create -n barcodemapper --override-channels -c conda-forge -c bioconda \
  python=3.12 snakemake=9.27.0 fastp bowtie2 samtools \
  pandas pyyaml biopython pysam plotly pip
conda activate barcodemapper
python -m pip install --no-deps \
  https://github.com/alberdilab/barcodemapper/releases/download/v1.2.2/barcodemapper-1.2.2-py3-none-any.whl
barcodemapper --version
```

## From source

The provided Conda environment pins all workflow tools and Python dependencies:

```sh
conda env create --file barcodemapper_env.yml
conda activate barcodemapper
python -m pip install --no-deps .
barcodemapper --version
```

The pinned environment is tested on Linux x86_64 (Mjolnir). fastp 1.4.0, Bowtie2 and samtools also publish native macOS builds, but other platforms are untested.

## Python-only install

`python -m pip install .` installs the Python dependencies. `snakemake`, `fastp`, `bowtie2`, and `samtools` must also be on `PATH`. Database builds and analyses without a current Bowtie2 index also need `bowtie2-build` (included with Bowtie2). Slurm runs additionally need the Snakemake Slurm executor plugin.

:::{note}
The short `pip install barcodemapper` and `conda install -c bioconda barcodemapper` commands require publication to PyPI and Bioconda, respectively, and are not yet available.
:::
