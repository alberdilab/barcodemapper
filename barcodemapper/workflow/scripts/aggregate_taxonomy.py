import pandas as pd
from collections import defaultdict

# Read input and output file paths from Snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Load taxonomic classifications (ignoring ReadID)
df = pd.read_csv(input_file, sep="\t")

# Count occurrences of each taxonomy, including missing assignments.
taxonomy_counts = df["Taxonomy"].fillna("").value_counts().to_dict()

# Initialize taxonomic hierarchy
taxonomy_tree = defaultdict(int)

def resolved_prefix(taxon):
    """Keep the highest contiguous ranks that have names."""
    levels = []
    for level in taxon.split(";"):
        if not level or level.endswith("__"):
            break
        levels.append(level)
    return ";".join(levels) if levels else "Unclassified"

# Process each taxon, ensuring parents accumulate child counts
for taxon, count in taxonomy_counts.items():
    taxon = resolved_prefix(taxon)

    # Assign count to the specific taxon
    taxonomy_tree[taxon] += count

    # Propagate counts to all parent taxonomic levels
    tax_levels = taxon.split(";")
    for i in range(1, len(tax_levels)):
        parent_taxon = ";".join(tax_levels[:i])
        taxonomy_tree[parent_taxon] += count

# Convert to DataFrame
agg_df = pd.DataFrame(list(taxonomy_tree.items()), columns=["Taxonomy", "Count"])

# Sort taxonomy alphabetically
agg_df = agg_df.sort_values(by="Taxonomy")

# Save output
agg_df.to_csv(output_file, sep="\t", index=False)
