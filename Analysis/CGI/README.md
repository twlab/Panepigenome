# CpG island annotation and methylation analysis

This directory contains code for identifying CpG islands (CGIs) from haplotype-resolved genome assemblies and for downstream analyses of CpG methylation within CGI-related genomic regions.

## Overview

CpG islands are identified from unmasked genome assembly sequences and converted to BED format for downstream pan-epigenomic analyses.

The workflows in this directory include:

- identification of CpG islands from haplotype-resolved genome assemblies;
- annotation of CGI shores and shelves; and
- calculation of mean and median CpG methylation stratified by genomic feature and variant class.

## CpG island annotation

CpG islands are identified from unmasked genome assembly sequences using `cpg_lh`.

Genome assemblies are first converted from FASTA to 2bit format using `faToTwoBit` and then converted back to unmasked FASTA sequences using `twoBitToFa -noMask` before CGI identification.

CpG islands identified by `cpg_lh` are converted to BED format for downstream analyses.

CGI shores are defined as regions extending up to 2 kb from CGI boundaries, excluding the CGIs themselves. CGI shelves are defined as regions 2–4 kb from CGI boundaries, excluding both CGIs and CGI shores.

Typical output files include:

```text
*.CGI.bed
*.CGIshores.bed
*.CGIshelves.bed

## Scripts

### CGI annotation workflow

The CGI annotation workflow identifies CpG islands from unmasked genome assembly sequences and generates genomic intervals for CpG islands, CGI shores and CGI shelves.

Major steps include:

1. conversion of FASTA sequences to 2bit format;
2. extraction of unmasked genome sequences;
3. identification of CpG islands using `cpg_lh`;
4. conversion of CGI annotations to BED format;
5. generation of CGI shores using regions extending 2 kb from CGI boundaries; and
6. generation of CGI shelves using regions extending 2–4 kb from CGI boundaries.

### `median_cpg_methylation.sh`

This script calculates mean and median CpG methylation for a specified CGI-related genomic feature and variant class.

The genomic feature and variant class are specified using:

```bash
fea="CGI"
fe="SV"
