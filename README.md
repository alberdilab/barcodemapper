# DietMapper

**DietMapper** is a software for quantifying plant, fungi and animal DNA by mapping shotgun sequencing reads against taxonomically annotated marker genes, including ITS, COI, 16S and others. It is useful for detecting and quantifying the relative proportion of dietary items from intestinal contents of faecal samples, providing the higher possible taxonomic resolution through resolving ambiguous mappings.

## 1. Installation

DietMapper can be installed with its python dependencies (but not snakemake, fastp, samtools and bowtie2) directly from this repository using **pip**.

```sh
pip install git+https://github.com/alberdilab/dietmapper.git
dietmapper -h
```

### Dependencies

#### Software

- snakemake=8.25.3
- fastp=0.23.1
- samtools=1.21
- bowtie2=2.4.2

#### Python modules

- numpy
- pandas
- argparse
- PyYAML
- requests
- plotly

If you don't have these dependencies installed and added to your path, you can create a custom ready-to-use environment for DietMapper using conda.

```sh
wget https://raw.githubusercontent.com/alberdilab/dietmapper/main/dietmapper_env.yml
conda env create --file dietmapper_env.yml
conda activate dietmapper
dietmapper -h
```

## 2. Source databases

There are two main modes for running DietMapper: **a)** using the default DietMapper database or **b)** generating a custom database from raw BOLD and UNITE databases.

### DietMapper database

This is the standard database containing marker genes sequences of plants, basidiomycota fungi, and animals. It was generated in **March 2025** with the most updated versions of BOLD (2025-03-14) and UNITE (*All eukaryotes; 2025-02-19) databases at the time. The database can be stored in any directory, and can be used in DietMapper by pointing towards the main fasta file `-d directory/to/database/dietmapper_db_202503.fa`.

```sh
cd dietmapper_db
wget https://sid.erda.dk/share_redirect/FGathhqKb5/dietmapper_db_202503.tar.gz
tar -xzvf dietmapper_db_202503.tar.gz
```

### Custom database

If you want to use different versions of the databases (e.g., a more recent update) or a different taxonomic breadth (e.g., exclude certain animals from the database), you can download the original databases and create a custom DietMapper database.

#### BOLD

1. Login to BOLDSystems (you need a user account).
2. Visit https://bench.boldsystems.org/index.php/datapackages/Latest
3. Obtain the link of the FASTA version of the database. Note that when making the download request a unique download link like the one shown below will be generated, so this step cannot be reproduced directly.
4. Download the FASTA file (gz compression) and decompress it.

```sh
cd dietmapper_db
curl -L -o BOLD_Public.14-Mar-2025.fa.gz "https://bench.boldsystems.org/index.php/API_Datapackage/fasta?id=BOLD_Public.14-Mar-2025&uid=167dcd55552bc4"
gunzip BOLD_Public.14-Mar-2025.fa.gz
```
#### UNITE

1. Visit: https://unite.ut.ee/repository.php
2. Obtain the link of the FASTA version of the All Eukaryotes database.
3. Download the FASTA file (tgz compression) and decompress it.

```sh
cd dietmapper_db
wget -O sh_general_release_dynamic_s_all_19.02.2025.tgz https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz
tar xvf sh_general_release_dynamic_s_all_19.02.2025.tgz
```

## 3. DietMapper usage

### With existing DietMapper database

Just provide the input data (-i), the path to the DietMapper database (-d) and the output (-o) file.

```sh
dietmapper -i path/to/reads_folder -d dietmapper_db_202503.fa -o final_file.txt
```

Or alternatively, add your forward (-1) and reverse (-2) sequencing reads as comma-separated lists.

```sh
dietmapper -1 sample1_1.fq.gz,sample2_1.fq.gz -2 sample1_2.fq.gz,sample2_2.fq.gz -d dietmapper_db_202503.fa -o final_file.txt
```

Note that while DietMapper can be run on individual samples, it is primarily designed for batch processing to generate combined output files (taxonomy file and sunburst plot).

### With original BOLD and UNITE databases

Provide the paths to the original BOLD (-b) and UNITE (-u) databases, as well as the destination of the DietMapper database (-d). Optionally, limit the taxa to be included in the DietMapper database using -x for limiting BOLD taxa and -y for limiting UNITE taxa. Once the DietMapper database has been generated, you will only need to use (-d).

```sh
dietmapper -i path/to/readsdir -b dietmapper_db/BOLD_Public.14-Mar-2025.fa -u dietmapper_db/sh_general_release_dynamic_s_all_19.02.2025.fasta -d dietmapper_db/dietmapper_db_202503.fa -x k__Animalia -y k__Viridiplantae,p__Basidiomycota -o dietmapper_results.txt
```

If you only want to create the database without running any analysis, use the --build argument without input and output arguments besides the databases.

```sh
dietmapper --build -b dietmapper_db/BOLD_Public.14-Mar-2025.fa -u dietmapper_db/sh_general_release_dynamic_s_all_19.02.2025.fasta -d dietmapper_db/dietmapper_db_202503.fa -x k__Animalia -y k__Viridiplantae,p__Basidiomycota
```

### On a computational cluster

When processing multiple samples on a computational cluster, it is recommended to use a screen session and SLURM to run the pipeline efficiently—optimising resource usage and adhering to the job queue system.

```sh
screen -S dietmapper
dietmapper -i path/to/readsdir -d dietmapper_db_202503.fa -o final_file.txt --slurm
```

## 4. Output files

### Taxonomy table

The main output of DietMapper is a quantitative taxonomy table that displays the number of reads mapped to each taxonomic level. Read counts are aggregated from lower to higher taxonomic ranks, providing the most specific classification possible—typically at the species level—when a read is uniquely mapped to a species-annotated reference. In cases where reads map to multiple taxa, the pipeline assigns the lowest common taxonomic level.

```
Taxonomy        SAMPLE1        SAMPLE2
k__Fungi        244     2
k__Fungi;p__Basidiomycota       244     2
k__Fungi;p__Basidiomycota;c__Tremellomycetes    223     0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Holtermanniales 1       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Holtermanniales;f__Holtermanniaceae     1       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Holtermanniales;f__Holtermanniaceae;g__Holtermannia     1       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Holtermanniales;f__Holtermanniaceae;g__Holtermannia;s__Holtermannia_saccardoi   1       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales     182     0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae  6       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Kwoniella     6       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Kwoniella;s__Kwoniella_shandongensis  6       0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae     59      0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae;g__Tremella 59      0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae;g__Tremella;s__Tremella_cheejenii   23      0
k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Tremellaceae;g__Tremella;s__Tremella_mesenterica 36      0
k__Viridiplantae        23      19
k__Viridiplantae;p__Anthophyta  23      12
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae       20      12
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales    19      9
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Betulaceae      0       8
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Betulaceae;g__Ostrya    0       7
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Betulaceae;g__Ostrya;s__Ostrya_carpinifolia     0       2
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae        13      0
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae;g__Castanea    4       0
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae;g__Castanea;s__Castanea_dentata        1       0
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae;g__Castanea;s__Castanea_seguinii       1       0
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae;g__Lithocarpus 1       0
k__Viridiplantae;p__Anthophyta;c__Eudicotyledonae;o__Fagales;f__Fagaceae;g__Lithocarpus;s__Lithocarpus_corneus  1       0
```
[Note that many lines have been removed from the below example]

### Sunburst plot

An optional second output of DietMapper is a HTML file containing sunburst plots of the quantitative taxonomic annotations.
