#!/bin/bash

id="HG00000"
sampleid="m00000"
dir="/storage1/fs1/ccdg2/Active/analysis/zdong/HPRCv2/QC/step2-aligment_QC/Align_asm/${id}"
mapq=10

# ----- Find maternal file -----
mat_candidate=$(ls "$dir/${id}_mat.${sampleid}"*.cram 2>/dev/null | head -1)

if [[ -n "$mat_candidate" && -e "$mat_candidate" ]]; then
    mat=$(basename "$mat_candidate")
else
    mat_candidate=$(ls "$dir/${id}_hap2.${sampleid}"*.cram 2>/dev/null | head -1)
    mat=$(basename "$mat_candidate")
fi


# ----- Find paternal file -----
pat_candidate=$(ls "$dir/${id}_pat.${sampleid}"*.cram 2>/dev/null | head -1)

if [[ -n "$pat_candidate" && -e "$pat_candidate" ]]; then
    pat=$(basename "$pat_candidate")
else
    pat_candidate=$(ls "$dir/${id}_hap1.${sampleid}"*.cram 2>/dev/null | head -1)
    pat=$(basename "$pat_candidate")
fi


hap1ref=$(ls /storage1/fs1/ccdg2/Active/analysis/HPRC_Release_2_Data/Release_2_Assemblies/${id}/${id}_hap1*.fa.gz 2>/dev/null | head -1)
if [[ -z "$hap1ref" ]]; then
    hap1ref=$(ls /storage1/fs1/ccdg2/Active/analysis/HPRC_Release_2_Data/Release_2_Assemblies/${id}/${id}_pat*.fa.gz 2>/dev/null | head -1)
fi


hap2ref=$(ls /storage1/fs1/ccdg2/Active/analysis/HPRC_Release_2_Data/Release_2_Assemblies/${id}/${id}_hap2*.fa.gz 2>/dev/null | head -1)
if [[ -z "$hap2ref" ]]; then
    hap2ref=$(ls /storage1/fs1/ccdg2/Active/analysis/HPRC_Release_2_Data/Release_2_Assemblies/${id}/${id}_mat*.fa.gz 2>/dev/null | head -1)
fi


### Extract required values

samtools view -@4 --input-fmt-option reference=$hap1ref "$dir/$pat" | awk 'BEGIN {
    OFS = "\t";
}
{
    seq = $9;
    mapq = $5;
    mg = rq = ".";

    for (i = 12; i <= NF; i++) {
        if ($i ~ /^mg:f:/) {
            split($i, a, ":"); mg = a[3];
        } else if ($i ~ /^rq:f:/) {
            split($i, a, ":"); rq = a[3];
        }
    }

    if (seq ~ /^[0-9]+$/ && seq > 0) {
        diff = length($10) - seq; 
        score = 1 - (diff < 0 ? -diff : diff) / length($10);
    } else {
        score = "NA";
    }

    print $1, mg, score, rq, mapq;
}' > "${pat}_original.log" &

samtools view -@4 --input-fmt-option reference=$hap2ref "$dir/$mat" | awk 'BEGIN {
    OFS = "\t";
}
{
    seq = $9;
    mapq = $5;
    mg = rq = ".";

    for (i = 12; i <= NF; i++) {
        if ($i ~ /^mg:f:/) {
            split($i, a, ":"); mg = a[3];
        } else if ($i ~ /^rq:f:/) {
            split($i, a, ":"); rq = a[3];
        }
    }

    if (seq ~ /^[0-9]+$/ && seq > 0) {
        diff = length($10) - seq; 
        score = 1 - (diff < 0 ? -diff : diff) / length($10);
    } else {
        score = "NA";
    }

    print $1, mg, score, rq, mapq;
}' > "${mat}_original.log" &
wait

cp "${pat}_original.log"  "${pat}_original_combined.log" &
cp "${mat}_original.log"  "${mat}_original_combined.log" & 
wait

### The base quality is filtered using the 'rq' tag, which may be applied during Circular Consensus Sequence (CCS) generation with the '--min-rq' option. However, we should verify that each file has applied this filtering properly.

if [[ -f "${pat}_original_combined.log" && -f "${mat}_original_combined.log" ]]; then
  awk -v mapq="$mapq" '$4 >= 0.9 && $5 >= mapq { print $1, $2, $3, $6 }' "${pat}_original_combined.log" > "$pat.log" &
  awk -v mapq="$mapq" '$4 >= 0.9 && $5 >= mapq { print $1, $2, $3, $6 }' "${mat}_original_combined.log" > "$mat.log" &
else
  awk -v mapq="$mapq" '$4 >= 0.9 && $5 >= mapq { print $1, $2, $3, 0 }' "${pat}_original.log" > "$pat.log" &
  awk -v mapq="$mapq" '$4 >= 0.9 && $5 >= mapq { print $1, $2, $3, 0 }' "${mat}_original.log" > "$mat.log" &
fi

wait

### Remove multi-mapped reads, keeping only the best-mapped one
sort -k1,1 -k2,2n -k3,3 "$pat.log" | awk '
{
    id = $1
    val2 = $2 + 0
    val3 = $3 + 0
    val4 = $4 + 0

    if (!(id in max2)) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    } else if (val2 * val3 > max2[id] * max3[id]) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    } else if (val2 * val3 == max2[id] * max3[id] && val2 > max2[id]) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    }
}
END {
    for (id in max2)
        print id, max2[id], max3[id], max4[id]
}
' | sort > "$pat.sort.log" &

sort -k1,1 -k2,2n -k3,3 "$mat.log" | awk '
{
    id = $1
    val2 = $2 + 0
    val3 = $3 + 0
    val4 = $4 + 0

    if (!(id in max2)) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    } else if (val2 * val3 > max2[id] * max3[id]) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    } else if (val2 * val3 == max2[id] * max3[id] && val2 > max2[id]) {
        max2[id] = val2
        max3[id] = val3
        max4[id] = val4
    }
}
END {
    for (id in max2)
        print id, max2[id], max3[id], max4[id]
}
' | sort > "$mat.sort.log" &
wait

### compare mg values in mapping to pat and mat
awk -v sampleid="${id}.$sampleid" '
FNR==NR {
    a2[$1] = $2
    a3[$1] = $3
    a4[$1] = $4
    next
}
{
    b2[$1] = $2
    b3[$1] = $3
    b4[$1] = $4
}
END {
    for (k in a2) {
        if (k in b2) {
            diff2 = a2[k] - b2[k]
            diff3 = a3[k] - b3[k]
            print k, a2[k], b2[k], diff2, a3[k], b3[k], diff3, a4[k], b4[k] >> (sampleid ".matched.log")
        }
    }

    for (k in a2) {
        if (!(k in b2)) {
            print k, a2[k], a3[k], a4[k] >> (sampleid ".unique_file1.log")
        }
    }

    for (k in b2) {
        if (!(k in a2)) {
            print k, b2[k], b3[k], b4[k] >> (sampleid ".unique_file2.log")
        }
    }
}
' "$pat.sort.log" "$mat.sort.log"


# Move only if files exist
[[ -f ${id}.${sampleid}.matched.log ]] && mv ${id}.${sampleid}.matched.log "$pat.$mat.matched.log"
[[ -f ${id}.${sampleid}.unique_file1.log ]] && mv ${id}.${sampleid}.unique_file1.log "$pat.unique.log"
[[ -f ${id}.${sampleid}.unique_file2.log ]] && mv ${id}.${sampleid}.unique_file2.log "$mat.unique.log"

### get reads that are uniquely mapped to hap1 and hap2
awk -v pat="$pat" -v mat="$mat" -v sampleid="${id}.$sampleid" '
{
    if ($2 * $5 > $3 * $6)
        print > (pat ".hap1_better.log")
    else if ($2 * $5 < $3 * $6)
        print > (mat ".hap2_better.log")
    else if ($2 * $5 == $3 * $6 && $4 > 0)
        print > (pat ".hap1_better.log")
    else if ($2 * $5 == $3 * $6 && $4 < 0)
        print > (mat ".hap2_better.log")
    else
        print > (sampleid ".both_better.log")
}' "$pat.$mat.matched.log"

shuf ${id}.${sampleid}.both_better.log | split -d -n 2 - ${id}.${sampleid}.part_

# Build file list for paternal
pat_files=()
[[ -f $pat.unique.log ]] && pat_files+=("$pat.unique.log")
[[ -f $pat.hap1_better.log ]] && pat_files+=("$pat.hap1_better.log")
[[ -f ${id}.${sampleid}.part_00 ]] && pat_files+=("${id}.${sampleid}.part_00")

cat "${pat_files[@]}" | awk '{print $1}' | sort > "$pat.better.txt"

# Build file list for maternal
mat_files=()
[[ -f $mat.unique.log ]] && mat_files+=("$mat.unique.log")
[[ -f $mat.hap2_better.log ]] && mat_files+=("$mat.hap2_better.log")
[[ -f ${id}.${sampleid}.part_01 ]] && mat_files+=("${id}.${sampleid}.part_01")

cat "${mat_files[@]}" | awk '{print $1}' | sort > "$mat.better.txt"

### extract sub-bam files
samtools view -@4 -N $pat.better.txt -o $pat.better.bam \
    --input-fmt-option reference=$hap1ref "$dir/$pat" &
samtools view -@4 -N $mat.better.txt -o $mat.better.bam \
    --input-fmt-option reference=$hap2ref "$dir/$mat" &
wait


