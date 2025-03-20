# dietscan
Quantification of dietary items from animal-derived metagenomes 

## Source databases

### BOLD

1. Login to BOLDSystems
2. Go to: https://bench.boldsystems.org/index.php/datapackages/Latest
3. Download data. Note that unique download codes will be generated, so this step cannot be reproduced directly.
```sh
curl -L -o BOLD.tar.gz "https://bench.boldsystems.org/index.php/API_Datapackage/fasta?id=BOLD_Public.14-Mar-2025&uid=167dba260d5969"
curl -L -o BOLD.fa.gz "https://bench.boldsystems.org/index.php/API_Datapackage/fasta?id=BOLD_Public.14-Mar-2025&uid=167dba260d5969"
```

Metazooa COI

```sh
wget https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz
```

### UNITE

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
```
