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
  * **CGI/** — CpG island annotation.
  * **CpG call/** — Phasing, alignment and methylation profiling.
  * **Population/** — Population-associated analyses.
  * **Genomic feature enrichment/** — Genomic feature enrichment analyses.
  * **Loci under directional selection/** — Analysis of loci under directional selection.
  * **mQTL/** — Methylation quantitative trait locus analyses.
  * **eQTM/** — Expression quantitative trait methylation analyses.
  * **Mediation and interaction/** — Mediation and interaction analyses.
  * **Clinical and pharmacogenomic annotation/** — Clinical and pharmacogenomic annotation.

* **Figures/**

  * Code used to reproduce the main, Extended Data and Supplementary figures.

## Software requirements

Analyses were performed using standard bioinformatics software and custom scripts written primarily in R, Python and Shell.

Major software used in the study includes:

* R
* minimap2
* pbmm2
* pb-CpG-tools
* CCS
* Primrose
* Jasmine
* bedtools

Software versions, command-line parameters and other dependencies are provided in the Methods of the manuscript and/or in the corresponding analysis scripts.

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
