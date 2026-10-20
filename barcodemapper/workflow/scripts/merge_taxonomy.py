import pandas as pd

# Get input and output files from Snakemake
input_files = list(snakemake.input)  # List of input taxonomic files
output_file = snakemake.output[0]  # Output file path
sample_names = list(snakemake.params.samples)  # Extract sample names from Snakemake wildcards

# Process each sample's aggregated taxonomy file
sample_frames = []
for file, sample_name in zip(input_files, sample_names):
    df = pd.read_csv(file, sep="\t")
    sample_frames.append(df.assign(_sample=sample_name))

# Combine sample counts into a wide table
if sample_frames:
    counts = pd.concat(sample_frames, ignore_index=True)
    counts = counts.drop_duplicates(subset=["Taxonomy", "_sample"], keep="last")
    merged_df = counts.pivot(index="Taxonomy", columns="_sample", values="Count")
else:
    merged_df = pd.DataFrame(index=pd.Index([], name="Taxonomy"))
merged_df = merged_df.reindex(columns=sample_names).infer_objects().fillna(0).astype(int)

# Reset index and rename columns
merged_df.index.name = "Taxonomy"
merged_df = merged_df.reset_index()

# Sort taxonomy alphabetically
merged_df = merged_df.sort_values(by="Taxonomy")

# Save final merged file
merged_df.to_csv(output_file, sep="\t", index=False)
