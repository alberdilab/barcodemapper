import pandas as pd
from collections import defaultdict
import sys

# Read input and output file paths from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Load TSV file
df = pd.read_csv(input_file, sep=" ", names=["ReadID", "Taxonomy", "Mismatches"])

# Group taxonomy assignments by read
read_taxonomy = defaultdict(list)

for _, row in df.iterrows():
    read_taxonomy[row["ReadID"]].append(row["Taxonomy"])

def find_finest_common_taxon(tax_list):
    """
    Given multiple taxonomic classifications, find the finest shared level.
    """
    split_taxonomies = [t.split(";") for t in tax_list]

    # Determine the common taxonomic levels
    min_length = min(map(len, split_taxonomies))  # Find shortest classification
    common_taxon = []

    for i in range(min_length):
        level = {tax[i] for tax in split_taxonomies}  # Collect unique values at this level
        if len(level) == 1:  # If all alignments share the same level, keep it
            common_taxon.append(list(level)[0])
        else:  # Stop at the last common level
            break

    return ";".join(common_taxon)

# Process each read
assigned_taxonomy = {read: find_finest_common_taxon(taxons) for read, taxons in read_taxonomy.items()}

# Convert to DataFrame
output_df = pd.DataFrame(assigned_taxonomy.items(), columns=["ReadID", "Taxonomy"])

# Save to TSV
output_df.to_csv(output_file, sep="\t", index=False)
