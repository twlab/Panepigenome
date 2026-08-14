#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO


def extract_cpg_to_bed(fasta_file, bed_file):
    cpg_count = 0

    with open(bed_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            chrom = record.id

            for i in range(len(seq) - 1):
                if seq[i:i + 2] == "CG":
                    # BED format: chrom, 0-based start, half-open end
                    out.write(f"{chrom}\t{i}\t{i + 2}\n")
                    cpg_count += 1

    return cpg_count


def main():
    parser = argparse.ArgumentParser(
        description="Identify all CpG dinucleotides in a FASTA file and output them in BED format."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("bed", help="Output BED file")
    args = parser.parse_args()

    try:
        cpg_count = extract_cpg_to_bed(args.fasta, args.bed)

        print(f"Output: {args.bed}")
        print(f"Total CpG sites: {cpg_count}")

    except FileNotFoundError:
        sys.exit(f"Error: input file not found: {args.fasta}")
    except Exception as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
