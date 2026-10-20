# BarcodeMapper

BarcodeMapper maps paired shotgun reads to taxonomically annotated COI and ITS reference sequences. It filters reads with fastp, aligns them with Bowtie2, and assigns each read ID to the lowest taxonomy shared by its accepted alignments. By default, one mapping pass produces four cumulative taxonomy count tables (ultrastrict, strict, standard and permissive), one self-contained HTML report, and one JSON run manifest.

:::{important}
Read counts describe marker-matching fragments. They are **not** estimates of organism abundance: marker copy number, reference coverage, read length, and mapping ambiguity affect the counts. See [](limitations.md).
:::

## Quick start

```sh
conda activate barcodemapper

# 1. Build the reference database (once)
barcodemapper --build --db-dir references   # creates references/bmdb_<YYYYMMDD>/

# 2. Analyse a directory of paired reads
barcodemapper -i reads --db-dir references -o results/taxonomy.tsv
```

The default run writes `results/taxonomy.ultrastrict.tsv`, `results/taxonomy.strict.tsv`, `results/taxonomy.standard.tsv`, `results/taxonomy.permissive.tsv`, `results/taxonomy.html` and `results/taxonomy.manifest.json`.

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
dependency-checks
changelog
```

BarcodeMapper source code is available under the [MIT license](https://github.com/alberdilab/barcodemapper/blob/main/LICENSE). Reference databases retain their providers' terms. Source: <https://github.com/alberdilab/barcodemapper>.
