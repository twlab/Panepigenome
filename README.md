# Human panepigenome

## About

This repository provides a collection of code for human panepigenomic analyses. The workflows include long-read DNA methylation quality control, projection of individual assembly coordinates to graph coordinates, graph-based interval-overlap analysis, and integrated representation of CpG copy number and methylation at var-CpGs. The repository also includes code used to reproduce the analyses and figures presented in the manuscript **“Human panepigenome represents epigenomic diversity.”**

## Repository structure

The repository is organized as follows:

* **QC/**

  * Long-read DNA methylation quality control.

* **CpG in graph coordinates/**

  * Projection of individual assembly coordinates to graph coordinates.

* **Intersection/**

  * Graph-based interval-overlap analysis.

* **Integrated representation/**

  * Integrated representation of CpG copy number and methylation at var-CpGs.

* **Analysis/**

  * **CpG count/** — Assembly-based CpG identification and counting.
  * **CGI/** — CpG island annotation and methylation analysis.
  * **CpG call/** — Haplotype-aware alignment and methylation profiling.
  * **Population/** — Population-associated analyses.
  * **Genomic feature enrichment/** — Genomic feature enrichment analyses.
  * **Loci under directional selection/** — Enrichment analysis at loci under directional selection.
  * **mQTL in var-CpG/** — Association between var-CpG methylation and genetic features.
  * **eQTM in var-CpG/** — Association between promoter var-CpG methylation and transcript expression.
  * **Mediation and interaction/** — Mediation and variant-by-methylation interaction analyses.
  * **Molecular QTLs/** — Enrichment of CpG-altering variants among molecular QTLs.
  * **Clinical and pharmacogenomic annotation/** — Clinical and pharmacogenomic annotation of CpG-altering variants.

* **Figures/**

  * Code used to reproduce the main, Extended Data and Supplementary figures.

## Software requirements

Analyses were performed using standard bioinformatics software and custom scripts written primarily in R, Python and Bash.

Major software used in the study includes:

| Software | Version | URL |
|---|---|---|
| bcftools | v1.15.1 | https://github.com/samtools/bcftools |
| bedops | v2.4.41 | https://github.com/bedops/bedops |
| BEDTools | v2.29.1 | https://bedtools.readthedocs.io/en/stable/ |
| CCS | v6.0.0, v6.2.0, v6.3.0, v6.4.0, v6.5.0, v7.0.0, v8.0.0, v8.0.1 | https://ccs.how/ |
| cpg_lh | — | http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/cpg_lh |
| DSS | v2.52.0 | https://www.bioconductor.org/packages/release/bioc/html/DSS.html |
| FastQTL | v2.0 | https://github.com/francois-a/fastqtl |
| LiftOver | — | https://genome.ucsc.edu/cgi-bin/hgLiftOver |
| MatrixEQTL | v2.3 | https://github.com/andreyshabalin/matrixeqtl |
| minimap2 | v2.26 | https://github.com/lh3/minimap2 |
| pb-CpG-tools | v2.3.2 | https://github.com/PacificBiosciences/pb-cpg-tools |
| pbmm2 | v1.13.1, v1.14.99 | https://github.com/PacificBiosciences/pbmm2 |
| PLINK | v1.90b7.4 | https://www.cog-genomics.org/plink/1.9/ |
| PLINK 2 | v2.00 | https://www.cog-genomics.org/plink/2.0/ |
| Primrose | v1.1.0, v1.2.0, v1.3.0, v1.4.0 | https://github.com/mattoslmp/primrose |
| Python 3 | v3.12.2, v3.12.10 | https://www.python.org/downloads/ |
| R | v4.4.1 | https://www.r-project.org/ |
| rust-fastQTL | — | https://github.com/huangnengCSU/rust-fastqtl |
| SAMtools | v1.21 | https://github.com/samtools/samtools |
| VCFtools | v0.1.16 | https://github.com/vcftools/vcftools |
| vg | v1.61.0 | https://github.com/vgteam/vg |

Where different software versions were used for data generated at different stages of the HPRC project, the relevant versions are specified in the corresponding scripts or workflow documentation.

## Data availability

The analyses use genomic, epigenomic and transcriptomic datasets described in the manuscript.

Primary HPRC datasets include:

- Haplotype-resolved genome assemblies and PacBio HiFi methylation data;
- Haplotype- and graph-based genetic variant calls;
- Iso-Seq data;
- Haplotype-resolved genome assembly annotations.

Accession numbers, URLs and additional information for the datasets are provided in the **Data availability** statement of the manuscript.

Large sequencing datasets and reference resources are not duplicated in this repository.

## Reproducibility

This repository contains the custom code and workflows used to reproduce the analyses and figures reported in the manuscript.

Scripts are organized according to the corresponding analyses and manuscript figures where applicable. Software versions, major parameters and statistical procedures are described in the Methods of the manuscript and/or accompanying scripts.

For analyses requiring high-performance computing, example job-submission scripts are provided where applicable.

## Figure reproduction

Code for generating the main, Extended Data and Supplementary figures is provided in the **Figures/** directory.

Where applicable, figure scripts are organized by figure or panel to facilitate reproduction of the results presented in the manuscript.

## Citation

If you use this code or analytical framework, please cite:

**Dong Z. et al. Human panepigenome represents epigenomic diversity.**

Journal and DOI information will be added upon publication.

## License

This repository is released under the **MIT License**. See the `LICENSE` file for details.

## Contact

For questions regarding the code or analyses, please contact:

**Zheng Dong**
Email: [dzheng.th@gmail.com](mailto:dzheng.th@gmail.com)
