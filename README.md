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

### UNITE

Plant and fungal ITS
