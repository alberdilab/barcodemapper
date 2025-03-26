import pandas as pd
from collections import defaultdict

# Get input and output files from Snakemake
input_files = list(snakemake.input)  # List of input taxonomic files
output_file = snakemake.output[0]  # Output file path
sample_names = list(snakemake.params.samples)  # Extract sample names from Snakemake wildcards

# Initialize dictionary for taxon counts
taxonomy_dict = defaultdict(lambda: defaultdict(int))

# Process each sample's aggregated taxonomy file
for file, sample_name in zip(input_files, sample_names):
    df = pd.read_csv(file, sep="\t")

    # Store counts for each taxon
    for _, row in df.iterrows():
        taxonomy = row["Taxonomy"]
        count = row["Count"]
        taxonomy_dict[taxonomy][sample_name] = count

# Convert dictionary to DataFrame
merged_df = pd.DataFrame.from_dict(taxonomy_dict, orient="index").fillna(0).astype(int)

# Reset index and rename columns
merged_df.index.name = "Taxonomy"
merged_df = merged_df.reset_index()

# Sort taxonomy alphabetically
merged_df = merged_df.sort_values(by="Taxonomy")

# Save final merged file
merged_df.to_csv(output_file, sep="\t", index=False)
