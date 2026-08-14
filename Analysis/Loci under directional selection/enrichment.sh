awk -F'\t' '
NR==FNR {
    key = "chr" $1 "\t" $2 "\t" $3 "\t" $4
    seen[key] = $0
    next
}
{
    key = $9 "\t" $10 "\t" $12 "\t" $13
    if (key in seen) {
        print seen[key] "\t" $0
    }
}
' Selection_Summary_Statistics_01OCT2025.hg38.txt /scratch/zdong/Projects/PanEpiG/V1-9/HPRC_VCF/DNAm/re.all.BOTH.txt > Both.overlap.txt

