### Assembly-based CpG identification and counting.

### `extract_CG_to_bed.py`

Identifies CpG dinucleotides from genome assembly sequences and converts their positions into BED format for downstream pan-epigenomic analyses.

#### Input
- Genome assembly sequence in FASTA format.

#### Output
- BED file containing the genomic coordinates of identified CpG sites.

CpG sites are identified by locating all occurrences of the dinucleotide `CG`, treating uppercase and lowercase sequences equivalently.

### Count CpG sites

The total number of identified CpG sites can be obtained by counting the number of lines in the output BED file:

```bash
wc -l example_output.bed
