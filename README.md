![BarcodeMapper logo](https://raw.githubusercontent.com/alberdilab/barcodemapper/main/images/barcodemapper.png)

# BarcodeMapper

[![Documentation](https://readthedocs.org/projects/barcodemapper/badge/?version=latest)](https://barcodemapper.readthedocs.io/en/latest/)
[![Package checks](https://github.com/alberdilab/barcodemapper/actions/workflows/package.yml/badge.svg)](https://github.com/alberdilab/barcodemapper/actions/workflows/package.yml)
[![Release](https://img.shields.io/github/v/release/alberdilab/barcodemapper)](https://github.com/alberdilab/barcodemapper/releases)
[![License: MIT](https://img.shields.io/github/license/alberdilab/barcodemapper)](LICENSE)
![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue)

BarcodeMapper maps paired shotgun reads to taxonomically annotated COI and ITS reference sequences. It filters reads with fastp, aligns them with Bowtie2, and assigns each read ID to the lowest taxonomy shared by its accepted alignments. It produces a cumulative taxonomy count table, a self-contained HTML report, and a JSON run manifest.

Read counts describe marker-matching fragments. They are **not** estimates of organism abundance. See [validation and limitations](https://barcodemapper.readthedocs.io/en/latest/limitations.html).

## Documentation

Full documentation is at **<https://barcodemapper.readthedocs.io>**:

- [Installation](https://barcodemapper.readthedocs.io/en/latest/installation.html)
- [Reference database](https://barcodemapper.readthedocs.io/en/latest/database.html)
- [Running an analysis](https://barcodemapper.readthedocs.io/en/latest/analysis.html)
- [Outputs](https://barcodemapper.readthedocs.io/en/latest/outputs.html)
- [Validation and limitations](https://barcodemapper.readthedocs.io/en/latest/limitations.html)

## Quick start

```sh
conda create -n barcodemapper --override-channels -c conda-forge -c bioconda \
  python=3.12 snakemake=8.27.1 fastp bowtie2 samtools \
  pandas pyyaml biopython pysam plotly pip
conda activate barcodemapper
python -m pip install --no-deps \
  https://github.com/alberdilab/barcodemapper/releases/download/v1.2.1/barcodemapper-1.2.1-py3-none-any.whl

# Build the reference database (once)
barcodemapper --build --db-dir references   # creates references/bmdb_<YYYYMMDD>/

# Analyse a directory of paired reads (sample_1.fq.gz / sample_2.fq.gz)
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv
```

This writes `results/taxonomy.tsv`, `results/taxonomy.html` and `results/taxonomy.manifest.json`.

## Citation and license

Cite the software release ([CITATION.cff](CITATION.cff)) and the reference providers (BOLD and UNITE). Changes are listed in [CHANGELOG.md](CHANGELOG.md).

BarcodeMapper source code is available under the [MIT license](LICENSE). Reference databases retain their providers' terms.
