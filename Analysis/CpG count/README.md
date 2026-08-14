### Assembly-based CpG identification and counting.

### `extract_CG_to_bed.py`

Identifies CpG dinucleotides from genome assembly sequences and converts their positions into BED format for downstream pan-epigenomic analyses.

#### Input
- Genome assembly sequence in FASTA format.

#### Output
- BED file containing the genomic coordinates of identified CpG sites.

CpG sites are identified by locating all occurrences of the dinucleotide `CG`, treating uppercase and lowercase sequences equivalently.
