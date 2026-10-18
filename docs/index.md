# BarcodeMapper

BarcodeMapper maps paired shotgun reads to taxonomically annotated COI and ITS reference sequences. It filters reads with fastp, aligns them with Bowtie2, and assigns each read ID to the lowest taxonomy shared by its accepted alignments. It produces a cumulative taxonomy count table, a self-contained HTML report, and a JSON run manifest.

:::{important}
Read counts describe marker-matching fragments. They are **not** estimates of organism abundance: marker copy number, reference coverage, read length, and mapping ambiguity affect the counts. See [](limitations.md).
:::

## Quick start

```sh
conda activate barcodemapper

# 1. Build the reference database (once)
barcodemapper --build -d references/barcodemapper.fa.gz

# 2. Analyse a directory of paired reads
barcodemapper -i reads -d references/barcodemapper.fa.gz -o results/taxonomy.tsv
```

```{toctree}
:maxdepth: 2
:caption: User guide

installation
database
analysis
outputs
limitations
```

```{toctree}
:maxdepth: 1
:caption: Project

development
changelog
```

BarcodeMapper source code is available under the [MIT license](https://github.com/alberdilab/barcodemapper/blob/main/LICENSE). Reference databases retain their providers' terms. Source: <https://github.com/alberdilab/barcodemapper>.
