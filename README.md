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

- fastp=0.23.1
- samtools=1.21
- bowtie2=2.4.2

If you don't have these dependencies installed, you can create a custom environment for DietScan using conda

```sh
wget 
conda env create --name dietscan --file dietscan_env.yml
```

## 2. Source databases

There are two main modes for running DietScan: using the default DietScan database or generating a custom database from raw BOLD and UNITE databases.

### DietScan database

This is the standard database containing ITS sequences of plants and basidiomycota fungi, and COI sequences of animals. It was generated in **March 2025** with the most updated versions of BOLD (2025-03-14) and UNITE (*All eukaryotes; 2025-02-19) databases at the time.

```sh
wget XXXXXXXXX
gunzip dietscan_db_202503.fa.gz
```

### Custom database

If you want to use other versions of the databases or a different taxonomic breadth, you can download the original databases and create a custom DietScan database.

#### BOLD

1. Login to BOLDSystems
2. Visit https://bench.boldsystems.org/index.php/datapackages/Latest
3. Download data. Note that unique download codes will be generated, so this step cannot be reproduced directly.

```sh
cd dietscan_db
curl -L -o BOLD_Public.14-Mar-2025.fa.gz "https://bench.boldsystems.org/index.php/API_Datapackage/fasta?id=BOLD_Public.14-Mar-2025&uid=167dcd55552bc4"
gunzip BOLD_Public.14-Mar-2025.fa.gz
```
#### UNITE

1. Visit: https://unite.ut.ee/repository.php
2. Download the General FASTA release for All Eukaryotes

```sh
cd dietscan_db
wget -O sh_general_release_dynamic_s_all_19.02.2025.tgz https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz
tar xvf sh_general_release_dynamic_s_all_19.02.2025.tgz
```

## 3. Analysis

### With existing DietScan database

Just provide the input data (-i), the path to the DietScan database (-d) and the output (-o) file.

```sh
dietscan -i path/to/readsdir -d dietscan_db_202305.fa -o final_file.txt
```

### With original BOLD and UNITE databases

Provide the paths to the original BOLD (-b) and UNITE (-u) databases, as well as the destination of the DietScan database (-d). Optionally, limit the taxa to be included in the DietScan database using -x for limiting BOLD taxa and -y for limiting UNITE taxa. Once the DietScan database has been generated, you will only need to use (-d)

```sh
dietscan -i path/to/readsdir -b dietscan_db/BOLD_Public.14-Mar-2025.fa -u dietscan_db/sh_general_release_dynamic_s_all_19.02.2025.fa -d dietscan_db/dietscan_db_202305.fa -x k__Animalia -y k__Viridiplantae,p__Basidiomycota -o DietScan_results.txt
```

### On a computational cluster

If using many samples in a computational cluster, it is better to use a screen session and slurm

```sh
screen -S dietscan
dietscan -i path/to/readsdir -d dietscan_db_202503.fa -o final_file.txt --slurm
```

### Using custom databases

```sh
wget -c custom_bold.fa [BOLD/URL]
wget -c custom_unite.fa [UNITE/URL]

screen -S dietscan
dietscan -i path/to/readsdir -b custom_bold.fa -u custom_unite.fa --bold_retain k__Animalia --unite_retain k__Viridiplantae,p__Basidiomycota -o final_file.txt --slurm
```

## Old scripts (to be removed when ready)

Plant and fungal ITS

## Database preparation

### Dereplicate BOLD database
Remove duplicated sequences

```sh
cat <<EOF > dereplicate_bold.sh
#!/bin/bash
#SBATCH --job-name=dereplicate_bold
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --time=6:00:00

python dereplicate_bold.py bold.fa bold_dep.fa
EOF

sbatch dereplicate_bold.sh
```

### Rename BOLD database
Unify header format and taxonomy

```sh
cat <<EOF > rename_bold.sh
#!/bin/bash
#SBATCH --job-name=rename_bold
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=6:00:00

python rename_bold.py bold_dep.fa bold_rename.fa
EOF

sbatch rename_bold.sh
```

### Rename UNITE database
Unify header format and taxonomy

```sh
cat <<EOF > rename_unite.sh
#!/bin/bash
#SBATCH --job-name=rename_unite
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=6:00:00

python rename_unite.py unite.fa unite_rename.fa
EOF

sbatch rename_unite.sh
```

### Merge both databases
With the option of retaining only desired taxa

```sh
cat <<EOF > merge_databases.sh
#!/bin/bash
#SBATCH --job-name=merge_databases
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=6:00:00

python merge_databases.py bold_rename.fa unite_rename.fa dietscan.fa --retain_bold k__Animalia --retain_unite k__Viridiplantae,k__Fungi
EOF

sbatch merge_databases.sh
```

### Combined
With the option of retaining only desired taxa

```sh
cat <<EOF > prepare_database.sh
#!/bin/bash
#SBATCH --job-name=merge_databases
#SBATCH --nodes=1
#SBATCH --partition=cpuqueue
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=2:00:00

python prepare_database.py -b bold.fa -u unite.fa -o dietscan.fa --retain_bold k__Animalia --retain_unite k__Viridiplantae,k__Fungi
EOF

sbatch prepare_database.sh
```
