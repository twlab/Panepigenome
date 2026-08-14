#!/bin/bash

i=HG00000

#tmprunfolder="/tmp/tmprunfolder_QJmCbO0uDNQV7DzQMYRw"
tmprunfolder=$(mktemp -d /tmp/tmprunfolder_XXXXXXXXXXXXXXXXXXXX)
#mkdir -p $tmprunfolder
cat /storage2/fs1/hprc/Active/zdong/Benchmark/RefGenome/Pangenome-GRCh38-full/hprc-v2.0-mc-grch38.full.gfa | grep "^W" | grep $i > $tmprunfolder/hprc-v2.0-mc-grch38.full.$i.W.gfa
cp /storage2/fs1/hprc/Active/zdong/Benchmark/RefGenome/Pangenome-GRCh38-full/hprc-v2.0-mc-grch38.full.S.gfa $tmprunfolder

hap=1
python3 lcpg2graph.py $tmprunfolder/hprc-v2.0-mc-grch38.full.S.gfa $tmprunfolder/hprc-v2.0-mc-grch38.full.$i.W.gfa /storage2/fs1/hprc/Active/zdong/V1-50/assembly-count/Bed/processed/processed_"$i".hap"$hap".bed $i $hap 20 > /storage2/fs1/hprc/Active/zdong/V1-50/assembly-count/Bed/Result/re.$i.hap$hap.txt 2> error.$i.hap$hap.err

hap=2
python3 lcpg2graph.py $tmprunfolder/hprc-v2.0-mc-grch38.full.S.gfa $tmprunfolder/hprc-v2.0-mc-grch38.full.$i.W.gfa /storage2/fs1/hprc/Active/zdong/V1-50/assembly-count/Bed/processed/processed_"$i".hap"$hap".bed $i $hap 20 > /storage2/fs1/hprc/Active/zdong/V1-50/assembly-count/Bed/Result/re.$i.hap$hap.txt 2> error.$i.hap$hap.err

rm -rf $tmprunfolder
