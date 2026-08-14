i=INS

python3 cpg_overlap_with_variant.py /storage2/fs1/hprc/Active/zdong/Benchmark/RefGenome/Pangenome-GRCh38-full/hprc-v2.0-mc-grch38.full.S.gfa variantSegments_"$i"3.txt final_counts_means_filtered.tsv  > re.all.$i.txt 2> error.$i.err
