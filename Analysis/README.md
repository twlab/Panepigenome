# Analysis

This directory contains code for downstream analyses of the human pan-epigenome dataset described in:

**Dong Z. et al. Human pan-epigenome represents epigenomic diversity.**

The analyses characterize CpG diversity, genomic context, population-associated variation, genetic regulation of DNA methylation, transcriptomic associations and potential clinical relevance of var-CpGs.

## Repository structure

- **CpG count/**
  - Assembly-based identification and quantification of CpG sites across haplotype-resolved genomes.

- **CGI/**
  - Annotation and analysis of CpG islands (CGIs), including CpG gain and loss within CGI regions.

- **CpG call/**
  - Processing of phased long-read alignments and DNA methylation calls.

- **Population/**
  - Population-associated analyses of CpG presence or absence and DNA methylation variation across continental groups and populations.

- **Genomic feature enrichment/**
  - Enrichment analyses of var-CpGs across genomic and regulatory features.

- **Loci under directional selection/**
  - Analysis of var-CpGs and associated loci showing evidence of population differentiation or directional selection.

- **mQTL/**
  - Identification and characterization of methylation quantitative trait loci (mQTLs) linking genetic variation to DNA methylation.

- **eQTM/**
  - Analysis of associations between DNA methylation and transcript expression.

- **Mediation and interaction/**
  - Mediation and genotype–methylation interaction analyses examining relationships among genetic variation, DNA methylation and transcript expression.

- **Clinical and pharmacogenomic annotation/**
  - Annotation of var-CpGs and associated variants using clinical and pharmacogenomic resources.

## Input data

These analyses use processed genomic, epigenomic and transcriptomic datasets generated or described in the manuscript, including:

- CpG presence or absence and methylation measurements;
- genetic variant genotypes;
- transcript expression data;
- population annotations;
- genomic and regulatory annotations; and
- external functional, clinical and pharmacogenomic resources.

Large primary datasets are not duplicated in this directory. Data access information is provided in the main repository README and the Data availability statement of the manuscript.

## Reproducibility

Scripts within each subdirectory correspond to the associated analysis described in the manuscript. Where applicable, analysis-specific README files provide information on input files, major parameters, software requirements and expected outputs.

Statistical methods, software versions and additional methodological details are described in the Methods of the manuscript.
