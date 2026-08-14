line1=***.fa
		/scratch/zdong/Softwares/cpglh/faToTwoBit $line1 $line1.2bit
		/scratch/zdong/Softwares/cpglh/twoBitToFa -noMask $line1.2bit stdout  \
                | /scratch/zdong/Softwares/cpglh/cpg_lh /dev/stdin 2> $line1.cpg_lh.err \
                |  awk '{$2 = $2 - 1; width = $3 - $2;  printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n", $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}' \
                | sort -k1,1 -k2,2n > $line1.CGI.bed

		samtools faidx $line1
		cut -f-2 $line1.fai > $line1.g
		bedtools slop -i $line1.CGI.bed -g $line1.g -b 2000 | bedtools merge \
		| bedtools subtract -a - -b $line1.CGI.bed > $line1.CGIshores.bed

		bedtools slop -i $line1.CGI.bed -g $line1.g -b 4000 | bedtools merge \
		| bedtools subtract -a - -b $line1.CGI.bed | bedtools subtract -a - -b $line1.CGIshores.bed > $line1.CGIshelves.bed
        #done < b.log


