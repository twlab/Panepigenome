bedtools intersect -a matched_meth.bed -b promoter | sort | uniq > matched_meth.promoter.bed
bedtools intersect -a matched_meth.bed -b /scratch/zdong/Projects/PanEpiG/References/CGI/GRCh38.cpg_lh.CpG_island.bed | sort | uniq > matched_meth.CGI.bed


grep -vFxf <(cat matched_meth.promoter.bed matched_meth.CGI.bed) matched_meth.bed > matched_meth.bed.filtered
