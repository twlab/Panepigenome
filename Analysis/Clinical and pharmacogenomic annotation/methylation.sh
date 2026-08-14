awk 'NR==FNR {
    if ($0 !~ /^#/) key[$1 FS $2 FS $4 FS $5] = 1
    next
}
{
    k = $9 FS $10 FS $12 FS $13
    if (k in key) print
}' ../0001.vcf /scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/DNAm/re.all.BOTH.txt > matched_meth.txt

# Step 1: Extract SNP lines into a temporary file
grep "SNP" /scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/DNAm/re.all.BOTH.txt > re.all.BOTH.SNP.txt

# Step 2: Match with another VCF
awk 'NR==FNR {
    if ($0 !~ /^#/) key[$1 FS $2 FS $4 FS $5] = 1
    next
}
{
    k = $9 FS $10 FS $12 FS $13
    if (k in key) print
}' ../0001.vcf re.all.BOTH.SNP.txt > matched_SNP_meth.txt
