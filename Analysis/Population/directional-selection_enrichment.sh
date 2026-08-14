# Population-associated analyses

This directory contains code for population-associated analyses of CpG abundance, DNA methylation, var-CpG copy number, integrated var-CpG signals and genetic–epigenetic differentiation across human populations.

## Overview

The workflows in this directory characterize population-associated variation in the human pan-epigenome using complementary CpG-level representations:

- methylation at shared CpGs;
- CpG copy number at var-CpGs;
- integrated var-CpG signals combining CpG copy number and methylation;
- population-specific CpG and methylation differences;
- population structure and variance partitioning;
- directional-selection enrichment; and
- genetic and epigenetic differentiation.

Unless otherwise stated, analyses are restricted to autosomal CpGs or var-CpGs.

## Diversity across methylomes

For each haplotype-resolved methylome, CpG abundance is summarized as the total number of detected CpGs and genome-wide methylation as the median methylation level across all CpGs.

Summary statistics are calculated by:

- population;
- continental group; and
- sex, where applicable.

Methylation levels are represented as proportions ranging from 0 to 1.

To quantify between-population variability within each continental group, population-level median CpG abundance or methylation is first calculated. The median absolute deviation (MAD) of these population-level values is then calculated within each continental group.

Differences between two groups are assessed using two-sided Wilcoxon rank-sum tests, whereas comparisons involving more than two populations or continental groups use Kruskal–Wallis tests.

Additional tests include:

- Pearson’s chi-squared tests for categorical distributions;
- Fisher’s exact tests for genomic-feature enrichment; and
- Benjamini–Hochberg correction for multiple testing.

Unless otherwise specified, FDR < 0.05 is considered significant.

## Variance partitioning and population structure

Population-associated epigenomic variation is quantified using three complementary CpG-level representations.

### Shared-CpG methylation

Shared CpGs are defined as autosomal reference CpGs present in more than 99% of analyzed haplotypes.

The response variable is DNA methylation level.

### Var-CpG copy number

Var-CpGs are restricted to autosomal sites with intermediate presence frequencies:

```text
5% < AF < 95%
