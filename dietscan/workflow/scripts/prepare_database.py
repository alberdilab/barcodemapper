import argparse
from Bio.Seq import Seq
from Bio import SeqIO
import random

### STEP 1: DEREPLICATE BOLD ###
def dereplicate_bold(input_fasta):
    """Removes redundancy:
    - For identical sequences: keeps the one with the most complete taxonomy.
    - For identical headers: keeps the one with the longest sequence.
    """
    from collections import defaultdict
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # First pass: deduplicate by sequence (keep the best taxonomy)
    seq_map = {}
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(record.seq)
        header = record.description
        taxonomy = header.split('|')[3] if len(header.split('|')) > 3 else ''

        if seq in seq_map:
            existing_record = seq_map[seq]
            existing_taxonomy = existing_record.description.split('|')[3] if len(existing_record.description.split('|')) > 3 else ''
            if len(taxonomy) > len(existing_taxonomy):
                seq_map[seq] = record
        else:
            seq_map[seq] = record

    # Second pass: deduplicate by header (keep the longest sequence)
    header_map = {}
    for record in seq_map.values():
        header = record.description
        seq = str(record.seq)
        if header in header_map:
            existing_seq = str(header_map[header].seq)
            if len(seq) > len(existing_seq):
                header_map[header] = record
        else:
            header_map[header] = record

    return list(header_map.values())

### STEP 2: RENAME BOLD HEADERS ###
def rename_bold(records):
    """Renames BOLD FASTA headers with taxonomic prefixes, filling empty levels with just the prefix."""
    renamed_records = []

    for record in records:
        parts = record.id.split('|')

        # Ensure taxonomy information is available
        taxonomy = parts[3] if len(parts) > 3 else ''
        taxon_levels = taxonomy.split(',') if taxonomy else []

        # Define taxonomic prefixes (excluding subfamily)
        taxon_prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

        # Remove subfamily (6th element) if it exists
        if len(taxon_levels) > 7:
            del taxon_levels[5]

        # Fill missing levels with empty strings
        while len(taxon_levels) < len(taxon_prefixes):
            taxon_levels.append('')

        # Format taxonomy
        formatted_taxonomy = []
        for prefix, taxon in zip(taxon_prefixes, taxon_levels):
            if taxon and taxon.lower() != 'none':
                formatted_taxonomy.append(f"{prefix}{taxon}")
            else:
                formatted_taxonomy.append(prefix)

        # Construct new header
        new_header = f"{parts[0]}|{parts[1]}|{';'.join(formatted_taxonomy)}"

        # Modify the record with the new header
        record.id = new_header
        record.description = ""
        renamed_records.append(record)

    return renamed_records

### STEP 3: RENAME UNITE HEADERS ###
def rename_unite(input_fasta):
    """Renames UNITE FASTA headers by extracting accession number and taxonomy."""
    renamed_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        parts = record.description.split('|')

        if len(parts) >= 5:
            accession = parts[1]  # UDB accession
            taxonomy = parts[-1]  # Last field is taxonomy
            new_header = f"{accession}|ITS|{taxonomy}"

            # Modify the record with new header
            record.id = new_header
            record.description = ""
            renamed_records.append(record)
        else:
            print(f"Skipping malformed record: {record.description}")

    return renamed_records

### STEP 4: FILTER & MERGE ###
def filter_fasta(records, retain_list):
    """Filters FASTA records, keeping only those that match specified taxa."""
    retained_records = []

    for record in records:
        parts = record.id.split('|')

        if len(parts) < 3:
            continue  # Skip malformed headers

        accession = parts[0]
        gene_region = parts[1]
        taxonomy_string = parts[2]

        taxonomy_levels = taxonomy_string.split(';')

        # Retain only if at least one specified taxon is found in the taxonomy part
        retain = any(taxon in taxonomy_levels for taxon in retain_list)

        if retain:
            # Reconstruct the header
            record.id = f"{accession}|{gene_region}|{taxonomy_string}"
            record.description = ""
            retained_records.append(record)

    return retained_records

def merge_fasta(bold_records, unite_records, output_fasta, retain_bold, retain_unite):
    """Merges the BOLD and UNITE databases after filtering."""
    filtered_bold = filter_fasta(bold_records, retain_bold)
    filtered_unite = filter_fasta(unite_records, retain_unite)

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(filtered_bold + filtered_unite, output_handle, "fasta")

### MAIN FUNCTION ###
def main():
    parser = argparse.ArgumentParser(description="Process, rename, and merge BOLD & UNITE FASTA databases with filtering.")
    parser.add_argument("-b","--bold_fasta", required=True, help="Input BOLD FASTA file")
    parser.add_argument("-u","--unite_fasta", required=True, help="Input UNITE FASTA file")
    parser.add_argument("-o","--output_fasta", required=True, help="Output merged FASTA file")
    parser.add_argument("-x","--retain_bold", required=False, default="", help="Comma-separated taxa to retain from BOLD (e.g. 'o__Coleoptera,o__Lepidoptera')")
    parser.add_argument("-y","--retain_unite", required=False, default="", help="Comma-separated taxa to retain from UNITE (e.g. 'k__Fungi,k__Viridioplantae')")

    args = parser.parse_args()

    # Convert comma-separated inputs into lists
    retain_bold = args.retain_bold.split(',') if args.retain_bold else []
    retain_unite = args.retain_unite.split(',') if args.retain_unite else []

    # Process BOLD dataset: Dereplicate & Rename
    print("Processing BOLD database...")
    bold_records = dereplicate_bold(args.bold_fasta)
    bold_records = rename_bold(bold_records)

    # Process UNITE dataset: Rename
    print("Processing UNITE database...")
    unite_records = rename_unite(args.unite_fasta)

    # Merge databases with filtering
    print("Merging and filtering databases...")
    merge_fasta(bold_records, unite_records, args.output_fasta, retain_bold, retain_unite)

    print("Processing complete! Output saved to:", args.output_fasta)

if __name__ == "__main__":
    main()
