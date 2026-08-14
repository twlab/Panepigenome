## Integrated representation of CpG copy number and methylation at var-CpGs

This workflow integrates CpG copy number and DNA methylation into a single quantitative representation for each var-CpG.

CpG copy number is encoded as one of three states (0, 1 or 2), and methylation is represented as a continuous value from 0 to 100%. Each observation is converted into a three-element vector, with the methylation value assigned to the element corresponding to the observed copy number state and the remaining elements set to −1.

A convolutional encoder–decoder network compresses this representation into a one-dimensional latent value. The latent value accurately captures both CpG copy number state (AUROC = 1.00) and methylation level (Pearson’s r = 0.99).

This latent value is used as the integrated representation of CpG copy number and methylation in downstream analyses.


## Contact

For questions about this tool, please contact:

**Saifur Rahman Jony**
Email: [saifur@wustl.edu](mailto:saifur@wustl.edu)
