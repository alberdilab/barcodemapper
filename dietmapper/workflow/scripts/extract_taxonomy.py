import pysam
import re
import argparse

def parse_cigar(cigartuples):
    """Calculate number of covered bases from CIGAR operations that consume both query and reference."""
    # CIGAR operation codes:
    # 0 = M, 7 = =, 8 = X → all consume query & reference
    return sum(length for op, length in cigartuples if op in {0, 7, 8})

def process_bam(input_bam, output_path, max_nm=2, min_covered=50):
    with pysam.AlignmentFile(input_bam, "rb") as bam, open(output_path, 'w') as out:
        for read in bam:
            if read.is_unmapped:
                continue

            try:
                nm = read.get_tag("NM")
            except KeyError:
                continue  # NM tag not present

            if nm > max_nm:
                continue

            covered_bases = parse_cigar(read.cigartuples)
            if covered_bases < min_covered:
                continue

            rname = bam.get_reference_name(read.reference_id)
            rnamesplit = rname.split("|")
            reference = rnamesplit[0] if len(taxon) > 0 else ''
            marker = rnamesplit[1] if len(taxon) > 1 else ''
            taxonomy = rnamesplit[2] if len(taxon) > 2 else ''

            out.write(f"{read.query_name}\t{reference}\t{marker}\t{taxonomy}\t{nm}\t{covered_bases}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter BAM reads by NM tag and covered bases.")
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output text file")
    parser.add_argument("--max_nm", type=int, default=2, help="Maximum allowed edit distance (NM tag)")
    parser.add_argument("--min_covered", type=int, default=50, help="Minimum number of covered bases from CIGAR")

    args = parser.parse_args()
    process_bam(args.input, args.output, max_nm=args.max_nm, min_covered=args.min_covered)
