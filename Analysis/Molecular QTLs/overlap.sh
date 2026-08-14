#!/bin/bash
set -euo pipefail

# ======= CONFIG =======
THREADS=8
REF_FASTA=/path/to/reference.fa   # only needed if you normalize alleles, can ignore
VCF_BIALLELIC=/scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/biallellic.vcf
AF_MIN=0.01
AF_MAX=0.99
OUTPUT_DIR=overlap_out

mkdir -p $OUTPUT_DIR

# ======= STEP 4: Find overlap using bcftools isec =======
echo "Step 4: Finding overlaps..."
bcftools isec --threads $THREADS --collapse none \
    -p $OUTPUT_DIR -n=2 \
    eQTL.vcf.gz \
    /scratch/zdong/Projects/PanEpiG/V1-9/ClinVar_common/biallellic_AF01_99.vcf.gz

echo "Pipeline complete. Overlap VCF is in $OUTPUT_DIR/0003.vcf"

bcftools view -H overlap_out/0000.vcf | wc -l
