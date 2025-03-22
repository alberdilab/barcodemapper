# DietScan

Quantification of dietary items from animal-derived metagenomes.

## 1. Installation

DietScan can be installed directly from this repository using **pip**.

```sh
pip install git+https://github.com/alberdilab/dietscan.git

#if you need to uninstall it:
pip uninstall dietscan -y
```

### Dependencies

#### Software

- fastp=0.23.1
- samtools=1.21
- bowtie2=2.4.2

#### Python modules

- numpy
- pandas
- argparse
- PyYAML
- requests
- snakemake
- plotly

If you don't have these dependencies installed and added to your path, you can create a custom environment for DietScan using conda.

```sh
wget https://raw.githubusercontent.com/alberdilab/dietscan/main/dietscan_env.yml
conda env create --name dietscan --file dietscan_env.yml
```

## 2. Source databases

There are two main modes for running DietScan: **a)** using the default DietScan database or **b)** generating a custom database from raw BOLD and UNITE databases.

### DietScan database

This is the standard database containing marker genes sequences of plants, basidiomycota fungi, and animals. It was generated in **March 2025** with the most updated versions of BOLD (2025-03-14) and UNITE (*All eukaryotes; 2025-02-19) databases at the time. The database can be stored in any directory, and can be used in DietScan by pointing towards the main fasta file `-d directory/to/database/dietscan_db_202503.fa`.

```sh
cd dietscan_db
wget https://sid.erda.dk/share_redirect/FGathhqKb5/dietscan_db_202503.tar.gz
tar -xzvf dietscan_db_202503.tar.gz
```

### Custom database

If you want to use different versions of the databases (e.g., a more recent update) or a different taxonomic breadth (e.g., exclude certain animals from the database), you can download the original databases and create a custom DietScan database.

#### BOLD

1. Login to BOLDSystems (you need a user account).
2. Visit https://bench.boldsystems.org/index.php/datapackages/Latest
3. Obtain the link of the FASTA version of the database. Note that when making the download request a unique download link like the one shown below will be generated, so this step cannot be reproduced directly.
4. Download the FASTA file (gz compression) and decompress it.

```sh
cd dietscan_db
curl -L -o BOLD_Public.14-Mar-2025.fa.gz "https://bench.boldsystems.org/index.php/API_Datapackage/fasta?id=BOLD_Public.14-Mar-2025&uid=167dcd55552bc4"
gunzip BOLD_Public.14-Mar-2025.fa.gz
```
#### UNITE

1. Visit: https://unite.ut.ee/repository.php
2. Obtain the link of the FASTA version of the All Eukaryotes database.
3. Download the FASTA file (tgz compression) and decompress it.

```sh
cd dietscan_db
wget -O sh_general_release_dynamic_s_all_19.02.2025.tgz https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz
tar xvf sh_general_release_dynamic_s_all_19.02.2025.tgz
```

## 3. DietScan usage

### With existing DietScan database

Just provide the input data (-i), the path to the DietScan database (-d) and the output (-o) file.

```sh
dietscan -i path/to/reads_folder -d dietscan_db_202503.fa -o final_file.txt
```

Or alternatively, add your forward (-1) and reverse (-2) sequencing reads as comma-separated lists.

```sh
dietscan -1 sample1_1.fq.gz,sample2_1.fq.gz -2 sample1_2.fq.gz,sample2_2.fq.gz -d dietscan_db_202503.fa -o final_file.txt
```

Note that while DietScan can be run on individual samples, it is primarily designed for batch processing to generate combined output files (taxonomy file and sunburst plot).

### With original BOLD and UNITE databases

Provide the paths to the original BOLD (-b) and UNITE (-u) databases, as well as the destination of the DietScan database (-d). Optionally, limit the taxa to be included in the DietScan database using -x for limiting BOLD taxa and -y for limiting UNITE taxa. Once the DietScan database has been generated, you will only need to use (-d).

```sh
dietscan -i path/to/readsdir -b dietscan_db/BOLD_Public.14-Mar-2025.fa -u dietscan_db/sh_general_release_dynamic_s_all_19.02.2025.fasta -d dietscan_db/dietscan_db_202305.fa -x k__Animalia -y k__Viridiplantae,p__Basidiomycota -o DietScan_results.txt
```

If you only want to create the database without running any analysis, use the --build argument without input and output arguments besides the databases.

```sh
dietscan --build -b dietscan_db/BOLD_Public.14-Mar-2025.fa -u dietscan_db/sh_general_release_dynamic_s_all_19.02.2025.fasta -d dietscan_db/dietscan_db_202305.fa -x k__Animalia -y k__Viridiplantae,p__Basidiomycota
```

### On a computational cluster

When processing multiple samples on a computational cluster, it is recommended to use a screen session and SLURM to run the pipeline efficiently—optimising resource usage and adhering to the job queue system.

```sh
screen -S dietscan
dietscan -i path/to/readsdir -d dietscan_db_202503.fa -o final_file.txt --slurm
```
