cut -f1 ../shared_var-cpg.txt | tail -n +2 | sort -u > shared.txt

awk 'NR==FNR {set[$1]; next} $1 in set' shared.txt \
/scratch/zdong/Projects/PanEpiG/V1-9/Geno-meth-exp/IndMeth/G-M_var-cpg/Matrixeqtls/cpg.coordinate.bed | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $1}' > shared.bed

awk 'NR==FNR {set[$1]; next} $1 in set' ../all_var-cpg_methbyvariants.txt \
/scratch/zdong/Projects/PanEpiG/V1-9/Geno-meth-exp/IndMeth/G-M_var-cpg/Matrixeqtls/cpg.coordinate.bed | awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $3, $4, $1}' > background.bed

# /scratch/zdong/Projects/PanEpiG/V1-9/ClinVar_common/overlap_Both/DNAm/promoter
bedtools intersect -a shared.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/full_stack_ChromHMM_annotations/wgEncodeBroadHmmGm12878HMM.bed -wa -wb > matched_shared.CREs.bed
bedtools intersect -a background.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/full_stack_ChromHMM_annotations/wgEncodeBroadHmmGm12878HMM.bed -wa -wb > matched_background.CREs.bed

bedtools intersect -a shared.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/exon -wa -wb > matched_shared.exon.bed
bedtools intersect -a background.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/exon -wa -wb > matched_background.exon.bed

bedtools intersect -a shared.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/intron -wa -wb > matched_shared.intron.bed
bedtools intersect -a background.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/intron -wa -wb > matched_background.intron.bed

bedtools intersect -a shared.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/Lymphoblastoid/merged.bed -wa -wb > matched_shared.SE.bed
bedtools intersect -a background.bed -b /scratch/zdong/Projects/PanEpiG/V1-9/ENCODE/Lymphoblastoid/merged.bed -wa -wb > matched_background.SE.bed

