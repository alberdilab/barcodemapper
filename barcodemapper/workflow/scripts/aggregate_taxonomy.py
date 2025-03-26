import pandas as pd
from collections import defaultdict

# Read input and output file paths from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Load taxonomic classifications (ignoring ReadID)
df = pd.read_csv(input_file, sep="\t")

# Count occurrences of each taxonomy
taxonomy_counts = df["Taxonomy"].value_counts().to_dict()

# Initialize taxonomic hierarchy
taxonomy_tree = defaultdict(int)

# Function to check for empty taxonomy levels
def has_empty_levels(taxon):
    return any(level.endswith("__") for level in taxon.split(";"))

# Process each taxon, ensuring parents accumulate child counts
for taxon, count in taxonomy_counts.items():
    if has_empty_levels(taxon):
        continue  # Skip taxa with empty levels

    # Assign count to the specific taxon
    taxonomy_tree[taxon] += count

    # Propagate counts to all parent taxonomic levels
    tax_levels = taxon.split(";")
    for i in range(1, len(tax_levels)):
        parent_taxon = ";".join(tax_levels[:i])
        if has_empty_levels(parent_taxon):
            continue  # Skip partial paths with empty levels
        taxonomy_tree[parent_taxon] += count

# Convert to DataFrame
agg_df = pd.DataFrame(list(taxonomy_tree.items()), columns=["Taxonomy", "Count"])

# Sort taxonomy alphabetically
agg_df = agg_df.sort_values(by="Taxonomy")

# Save output
agg_df.to_csv(output_file, sep="\t", index=False)
