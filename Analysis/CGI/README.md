# CpG island annotation

This directory contains code for identifying and annotating CpG islands (CGIs) in haplotype-resolved genome assemblies and for downstream analyses of CpG gain, loss and methylation within CGI regions.

## Overview

CpG islands are annotated from genome assembly sequences and used to characterize the genomic context of CpGs and var-CpGs.

The workflows in this directory include:

- identification of CpG islands from haplotype-resolved assemblies;
- conversion of CGI annotations to BED format;
- overlap analysis between CGIs and CpG or var-CpG intervals;
- comparison of CpG gain and loss within CGI regions; and
- downstream methylation and enrichment analyses involving CGIs.

## Input

Typical inputs include:

- haplotype-resolved genome assembly sequences in FASTA format;
- CpG coordinates in BED format;
- var-CpG coordinates;
- genetic variant annotations; and
- CGI annotation files generated for each assembly.

## Output

Depending on the analysis, outputs may include:

- CGI coordinates in BED format;
- CpGs or var-CpGs overlapping CGI regions;
- counts of CpG gain and loss within CGIs;
- methylation measurements for CGI-associated CpGs; and
- summary tables used for downstream statistical analyses and figure generation.

## Genomic overlap

CGI and CpG intervals are intersected using genomic coordinates to identify CpGs located within CGI regions.

Where appropriate, complete overlap criteria are used to ensure that the CpG interval is contained within the annotated CGI.

## Downstream analyses

CGI annotations are used to evaluate:

- the distribution of CpG gain and loss across CGI regions;
- methylation patterns of CGI-associated var-CpGs;
- enrichment of var-CpGs within CGIs; and
- CGI-level effects associated with genetic variation.

## Reproducibility

Scripts in this directory correspond to the CGI-related analyses described in the manuscript.

Software versions, parameters and statistical procedures are provided in the corresponding scripts and in the Methods of the manuscript.
