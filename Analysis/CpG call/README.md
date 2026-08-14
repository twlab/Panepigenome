````markdown
# Phasing, alignment and methylation profiling

This directory contains code for haplotype-aware alignment of PacBio HiFi reads, read assignment to haplotype-resolved assemblies and single-base DNA methylation profiling.

## Overview

PacBio HiFi reads are aligned independently to the two haplotype-resolved assemblies available for each individual. Reads are then assigned to the haplotype providing stronger alignment support, followed by haplotype-specific DNA methylation calling.

For trio-phased assemblies, the two haplotypes correspond to paternal and maternal assemblies. For other samples, they are designated haplotypes 1 and 2.

The workflow includes:

- alignment of HiFi reads to both haplotype assemblies;
- extraction of alignment and read-quality metrics;
- filtering of low-quality reads and alignments;
- haplotype assignment based on relative alignment support;
- generation and sorting of haplotype-specific BAM files;
- single-base 5mC profiling using `pb-CpG-tools`; and
- filtering of methylation calls according to CpG presence and sequencing depth.

## Read alignment

PacBio HiFi reads are aligned independently to both haplotype-resolved genome assemblies using `minimap2`.

For each read alignment, the following information is extracted:

- read identifier;
- mapping quality;
- alignment score (`mg`);
- CCS read quality (`rq`);
- aligned read length; and
- total read length.

Reads are retained if:

```text
rq >= 0.9
mapping quality >= 10
````

For reads with multiple alignments within the same haplotype assembly, the best-supported alignment is retained.

## Haplotype assignment

Alignment support is calculated using the alignment score and read-length consistency:

```text
alignment support = mg × (aligned read length / total read length)
```

Reads aligning to both haplotypes are assigned to the haplotype with the stronger alignment support.

Reads with equivalent support between the two haplotypes are randomly assigned to one haplotype to maintain balanced coverage for downstream analyses.

The resulting reads are written to haplotype-specific BAM files.

## BAM processing

Haplotype-assigned BAM files are sorted using `SAMtools`.

Example:

```bash
samtools sort -o haplotype1.sorted.bam haplotype1.bam
samtools sort -o haplotype2.sorted.bam haplotype2.bam
```

BAM files can subsequently be indexed using:

```bash
samtools index haplotype1.sorted.bam
samtools index haplotype2.sorted.bam
```

## DNA methylation profiling

DNA methylation levels are quantified at single-base resolution using `pb-CpG-tools`.

The methylation pileup is generated using:

```bash
pb-CpG-tools aligned_bam_to_cpg_scores \
    --min-mapq 10 \
    --pileup-mode model \
    ...
```

The `model` pileup mode estimates CpG methylation levels from the `MM` and `ML` tags in PacBio HiFi reads.

Methylation profiling is performed separately for each haplotype-assigned BAM file.

## CpG and coverage filtering

To reduce bias caused by sequence differences between reads and haplotype assemblies, methylation values are retained only for CpG sites that are present in the corresponding haplotype assembly.

CpG methylation calls are additionally required to have a minimum sequencing depth of:

```text
10×
```

Only CpGs satisfying both criteria are retained for downstream pan-epigenomic analyses.

## Input

Typical inputs include:

* PacBio HiFi BAM or read files;
* two haplotype-resolved genome assemblies for each individual;
* CpG coordinates identified from each haplotype assembly; and
* sample and haplotype metadata.

## Output

The workflow generates:

* read alignments against each haplotype assembly;
* haplotype-assigned BAM files;
* sorted and indexed haplotype-specific BAM files;
* single-base CpG methylation calls;
* sequencing-depth information; and
* filtered haplotype-specific methylation profiles used in downstream analyses.

## Software requirements

Major software used in this workflow includes:

* `minimap2`
* `SAMtools` v1.21
* `pb-CpG-tools` v2.3.2
* Python and/or standard Unix command-line utilities used for read assignment and filtering

Exact software versions and parameters are provided in the corresponding scripts and in the Methods of the manuscript.

## Reproducibility

Scripts in this directory correspond to the phasing, haplotype-aware alignment and DNA methylation profiling procedures described in the manuscript.

The workflow retains the relationship between each HiFi read, its assigned haplotype and the corresponding assembly-specific CpG methylation measurements for downstream pan-epigenomic analyses.

