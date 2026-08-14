#!/bin/bash

# Define the population list
pops=(afr amr eas eur sas)

vcf="/scratch/zdong/Projects/PanEpiG/V1-9/Model_Ind/Methmatch/Var/merged_0label_filteredsex.vcf"
pop_dir="/scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/FST"

for pop in "${pops[@]}"; do
    echo "Calculating FST for ${pop} vs all other populations ..."

    # Pop1 = current population
    pop1="${pop_dir}/${pop}.txt"

    # Pop2 = all other populations combined
    other_pops=()
    for other in "${pops[@]}"; do
        if [[ "$other" != "$pop" ]]; then
            other_pops+=("${pop_dir}/${other}.txt")
        fi
    done
    pop2="${pop}_other.txt"
    cat "${other_pops[@]}" > "$pop2"

    # Run vcftools FST
    vcftools \
      --vcf "$vcf" \
      --not-chr chrX \
      --not-chr chrY \
      --not-chr chrM \
      --weir-fst-pop "$pop1" \
      --weir-fst-pop "$pop2" \
      --maf 0.05 \
      --out "${pop}_vs_rest"
done
