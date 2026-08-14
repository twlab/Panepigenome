````markdown
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
````

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

An example workflow is:

```bash
faToTwoBit input.fa input.fa.2bit

twoBitToFa -noMask input.fa.2bit stdout \
    | cpg_lh /dev/stdin 2> input.fa.cpg_lh.err \
    | awk '{
        $2 = $2 - 1;
        width = $3 - $2;
        printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n",
               $1, $2, $3, $5, $6, width, $6,
               width*$7*0.01, 100.0*2*$6/width, $7, $9);
      }' \
    | sort -k1,1 -k2,2n > input.fa.CGI.bed

samtools faidx input.fa
cut -f1-2 input.fa.fai > input.fa.genome

bedtools slop \
    -i input.fa.CGI.bed \
    -g input.fa.genome \
    -b 2000 \
    | bedtools merge \
    | bedtools subtract \
        -a - \
        -b input.fa.CGI.bed \
    > input.fa.CGIshores.bed

bedtools slop \
    -i input.fa.CGI.bed \
    -g input.fa.genome \
    -b 4000 \
    | bedtools merge \
    | bedtools subtract \
        -a - \
        -b input.fa.CGI.bed \
    | bedtools subtract \
        -a - \
        -b input.fa.CGIshores.bed \
    > input.fa.CGIshelves.bed
```

### `median_cpg_methylation.sh`

This script calculates mean and median CpG methylation for a specified CGI-related genomic feature and variant class.

The genomic feature and variant class are specified using:

```bash
fea="CGI"
fe="SV"
```

The script:

* reads sample and haplotype identifiers from `input.txt`;
* processes the corresponding methylation BED files;
* filters records by the selected variant class;
* removes duplicate CpG records using the first four columns as a unique identifier;
* calculates the mean CpG methylation for each sample-haplotype pair; and
* calculates the median CpG methylation for each sample-haplotype pair.

The script is:

```bash
fea="CGI"
fe="SV"
outfile="summary.${fe}.${fea}.mean.tsv"

echo -e "ID\tmean\tmedian" > "$outfile"

while read -r id hap; do
    [[ -z "$id" || "$id" =~ ^# ]] && continue

    over="${id}.${hap}.del.$fea.meth.bed"

    [[ ! -s "$over" ]] && {
        echo -e "${id}.${hap}\t0\t0" >> "$outfile"
        continue
    }

    read over_mean over_median < <(
        awk -F'\t' -v fe="$fe" '
            $18 == fe {
                key = $1 FS $2 FS $3 FS $4
                if (!seen[key]++) {
                    vals[++n] = $4
                    sum += $4
                }
            }
            END {
                if (n == 0) {
                    print 0, 0
                    exit
                }

                asort(vals)
                mean = sum / n

                if (n % 2)
                    median = vals[(n + 1) / 2]
                else
                    median = (vals[n / 2] + vals[n / 2 + 1]) / 2

                print mean, median
            }
        ' "$over"
    )

    echo -e "${id}.${hap}\t${over_mean}\t${over_median}" >> "$outfile"

done < input.txt
```

If no matching CpGs are identified for a sample-haplotype pair, both the mean and median are reported as `0`.

The output file is:

```text
summary.${fe}.${fea}.mean.tsv
```

with the following columns:

```text
ID    mean    median
```

## Input

Typical inputs include:

* haplotype-resolved genome assembly sequences in FASTA format;
* CGI, CGI shore and CGI shelf annotations in BED format;
* CpG methylation BED files;
* variant-class annotations; and
* sample and haplotype identifiers provided in `input.txt`.

An example `input.txt` is:

```text
sample1    1
sample1    2
sample2    1
sample2    2
```

where the first column represents the sample identifier and the second column represents the haplotype.

## Output

Outputs generated by these workflows include:

* CpG island coordinates;
* CGI shore coordinates;
* CGI shelf coordinates;
* mean CpG methylation for each sample-haplotype pair; and
* median CpG methylation for each sample-haplotype pair.

## Genomic features

The methylation-summary workflow can be applied to different CGI-related genomic features by changing the `fea` variable.

For example:

```bash
fea="CGI"
```

or:

```bash
fea="CpG_shore"
```

The corresponding input filenames are constructed from the sample identifier, haplotype and selected genomic feature.

## Variant classes

CpG methylation summaries can be stratified by genetic variant class using the `fe` variable.

For example:

```bash
fe="SNV"
```

or:

```bash
fe="SV"
```

Only records matching the selected variant class in column 18 of the input BED file are included in the methylation summary.

## Requirements

Major software and command-line utilities used in these workflows include:

* `faToTwoBit`
* `twoBitToFa`
* `cpg_lh`
* `samtools`
* `bedtools`
* `awk`
* `sort`
* standard Unix command-line utilities

## Reproducibility

Scripts in this directory correspond to CGI-related analyses described in the manuscript.

Software versions, parameters and additional methodological details are provided in the corresponding scripts and in the Methods of the manuscript.

```
