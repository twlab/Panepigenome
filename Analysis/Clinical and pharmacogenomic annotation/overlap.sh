#!/bin/bash
set -euo pipefail

# ======= CONFIG =======
THREADS=8

bcftools isec --threads $THREADS --collapse none \
    -p ./ -n=2 \
    ../clinvar_chr_fixed.vcf.gz \
    /scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/DNAm/BOTH_positions.vcf.sorted.gz
