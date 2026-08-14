#!/usr/bin/env python3
"""
extract_CG_to_bed.py

Identify all CpG dinucleotides in a FASTA file and output their
coordinates in BED format.

Usage:
    python3 extract_CG_to_bed.py input.fa output.bed
"""

import sys
import argparse
from Bio import SeqIO


def extract_cpg_to_bed(fasta_file, bed_file):
    with open(bed_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            chrom = record.id

            for i in range(len(seq) - 1):
                if seq[i:i + 2] == "CG":
                    # BED format: chrom, 0-based start, half-open end
                    out.write(f"{chrom}\t{i}\t{i + 2}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Identify all CpG dinucleotides in a FASTA file "
            "and output their coordinates in BED format."
        )
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("bed", help="Output BED file")
    args = parser.parse_args()

    try:
        extract_cpg_to_bed(args.fasta, args.bed)
        print(f"CpG extraction completed successfully. Output: {args.bed}")
    except FileNotFoundError:
        sys.exit(f"Error: input file not found: {args.fasta}")
    except Exception as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
