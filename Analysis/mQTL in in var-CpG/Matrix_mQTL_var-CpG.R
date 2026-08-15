# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

library(MatrixEQTL)

## Location of the package with the data files.
base.dir = '/scratch/zdong/Projects/PanEpiG/V1-9/Geno-meth-exp/IndMeth/G-M_var-cpg/Matrixeqtls';

## Settings

# Linear model to use
useModel = modelLINEAR;

# Genotype file name
SNP_file_name = paste(base.dir, "/var_aligned.bed", sep="");
snps_location_file_name = paste(base.dir, "/var.coordinate.bed", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/meth_aligned.bed", sep="");
gene_location_file_name = paste(base.dir, "/cpg.coordinate.bed", sep="");

# Covariates file name
covariates_file_name = paste(base.dir, "/cov_aligned.txt", sep="");

# Output file name
output_file_name_cis = "result_cis.txt";
output_file_name_tra = NULL;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 1e4;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      
snps$fileOmitCharacters = "NA"; 
snps$fileSkipRows = 1;          
snps$fileSkipColumns = 1;     
snps$fileSliceSize = 2000;    
snps$LoadFile(SNP_file_name);

snps$columnNames

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      
gene$fileOmitCharacters = "NA"; 
gene$fileSkipRows = 1;         
gene$fileSkipColumns = 1;       
gene$fileSliceSize = 2000;    
gene$LoadFile(expression_file_name);

gene$columnNames

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      
cvrt$fileOmitCharacters = "NA"; 
cvrt$fileSkipRows = 1;          
cvrt$fileSkipColumns = 1;      
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

## Plot the Q-Q plot of local and distant p-values

me$cis$ntests
plot(me)

