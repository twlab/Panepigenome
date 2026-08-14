### Set working directory
setwd("./")

### Load required packages
library(DSS)     
library(bsseq)   

### Define the data path
path <- file.path("./")

### Read files
dat1.1 = read.table(file.path(path, "example1_input_5.1and5.2_dmr.log"), header=T)
dat1.2 = read.table(file.path(path, "example2_input_5.1and5.2_dmr.log"), header=T)

### Select CpGs with coverage ≥ 10
dat1.1<-dat1.1[dat1.1$total_coverage>9,]
dat1.2<-dat1.2[dat1.2$total_coverage>9,]

### Columns were renamed according to the DSS package documentation
names(dat1.1) <- c("chr", "pos", "N", "X")
names(dat1.2) <- c("chr", "pos", "N", "X")

### Create an object of BSseq class, which is defined in bsseq Bioconductor package
BSobj = makeBSseqData(list(dat1.1, dat1.2),
                      c("C1","N1") )

### Smooth the methylation levels
dmlTest = DMLtest(BSobj, group1=c("C1"), group2=c("N1"), 
                  smoothing=TRUE,ncores=4)

### to detect CpGs with difference greater than 0.2
dmls2 = callDML(dmlTest, delta=0.2, p.threshold=0.001) # Delta represents the absolute value of the methylation change
head(dmls2)
write.table(dmls2, file = "DML_results_cov10_0.2_0.001.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

### Detect DMRs: delta > 0.2, min length = 50 bp, min CpGs = 3
dmls3 = callDMR(dmlTest, delta=0.2, p.threshold=0.01, minlen=50, minCG=3, dis.merge=50, pct.sig=0.5) # Parameter details are available at https://rdrr.io/bioc/DSS/man/callDMR.html
head(dmls3)
write.table(dmls3, file = "DMR_results_cov10_0.2_0.01.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

