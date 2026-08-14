LSF_DOCKER_PRESERVE_ENVIRONMENT=false
LSF_DOCKER_VOLUMES="/storage2/fs1/hprc/Active:/storage2/fs1/hprc/Active:/storage2/fs1/hprc/Active:/storage2/fs1/hprc/Active /storage1/fs1/ccdg2/Active:/storage1/fs1/ccdg2/Active"

i=name
dir=XXXXX (set this dirctory to your path that do the analysis)

### step 1: Alignment
bsub -oo $i.hg38.out -J pb$i.hg38 -G compute-hprc -q general -M 40G -n 2 -u zdong@wustl.edu -N  -R "select[mem>200000] rusage[mem=200000] span[hosts=1]" -a 'docker(shihuadong/biotools)' pbmm2 align /storage1/fs1/hprc/Active/zdong/Benchmark/RefGenome/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi $i.fofn $dir/$i.hg38.bam --preset HIFI --unmapped --num-threads 2 --log-level INFO --log-file $dir/pbmm2.$i.hg38.log --sample $i"

### step 2: Sort
bsub -oo $i.sort.out -J sort$i.hg38 -G compute-hprc -q general -M 40G -n 4 -u zdong@wustl.edu -N  -R "select[mem>40000 && tmp>5] rusage[mem=40000, tmp=5] span[hosts=1]" -a 'docker(shihuadong/biotools)' samtools sort -@ 4 $dir/$i.hg38.bam -o $dir/$i.hg38.sort.bam

### step 3: Index
bsub -oo $i.sort.out -J sort$i.hg38 -G compute-hprc -q general -M 40G -n 4 -u zdong@wustl.edu -N  -R "select[mem>40000 && tmp>5] rusage[mem=40000, tmp=5] span[hosts=1]" -a 'docker(shihuadong/biotools)' samtools index $dir/$i.hg38.sort.bam -@ 4

### step 4: CpG call
bsub -oo $i.cpgcall.out -J CpGcall$i -G compute-hprc -q general -M 40G -n 4 -u zdong@wustl.edu -N  -R "select[mem>40000] rusage[mem=40000,tmp=5] span[hosts=1]" -a 'docker(shihuadong/biotools)' /storage2/fs1/hprc/Active/zdong/Software/bin/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
  --bam $dir/$i.hg38.sort.bam \
  --output-prefix $dir/CpGcall/$i.model.pbmm2 \
  --model /storage2/fs1/hprc/Active/zdong/Software/bin/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite \
  --threads 4 \
  --min-mapq 10