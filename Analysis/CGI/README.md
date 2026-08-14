# CpG island annotation and methylation analysis

This directory contains code for identifying CpG islands (CGIs) from haplotype-resolved genome assemblies and for summarizing CpG methylation within CGI-related genomic regions.

## Overview

CpG islands are identified from unmasked genome assembly sequences and converted to BED format for downstream pan-epigenomic analyses.

The workflows in this directory include:

- identification of CpG islands from haplotype-resolved assemblies;
- annotation of CGI shores and shelves;
- overlap-based extraction of CpGs associated with CGI-related regions; and
- calculation of mean and median CpG methylation stratified by variant class.

## CpG island annotation

CpG islands are identified from unmasked genome assembly sequences using `cpg_lh`.

Genome assemblies are converted from FASTA to 2bit format and subsequently converted back to unmasked FASTA sequences using `twoBitToFa -noMask` before CGI identification.

CGI shores are defined as regions extending up to 2 kb from CGI boundaries, excluding the CGIs themselves. CGI shelves are defined as regions 2–4 kb from CGI boundaries, excluding both CGIs and CGI shores.

## Scripts

### CGI annotation script

Identifies CpG islands from unmasked genome assemblies and generates BED files for:

- CpG islands;
- CGI shores; and
- CGI shelves.

Typical output files include:

```text
*.CGI.bed
*.CGIshores.bed
*.CGIshelves.bed
