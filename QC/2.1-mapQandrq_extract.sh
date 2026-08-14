#!/bin/bash

dirr="/storage1/fs1/ccdg2/Active/analysis/zdong/HPRCv2/QC/step2-aligment_QC/Aligh_t2t"
input_file="example_input_2.1-mapQandrq_extract.log"
output_count="mapq_count.txt"

# Optional: header
echo -e "Sample\tTotal\t<60\t<50\t<40\t<30\trq<0.99" > "$output_count"

while read -r name; do
    file="${dirr}/${name}"
    
    samtools view -@2 "$file" | awk '
    BEGIN {
        OFS = "\t";
        prev = "";
        total = lt60 = lt50 = lt40 = lt30 = rq_below = 0;
    }
    {
        seq = $9 + 0;
        mapq = $5;
        mg = qs = qe = rq = ".";
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^mg:f:/) { split($i, a, ":"); mg = a[3]; }
            else if ($i ~ /^qs:i:/) { split($i, a, ":"); qs = a[3]; }
            else if ($i ~ /^qe:i:/) { split($i, a, ":"); qe = a[3]; }
            else if ($i ~ /^rq:f:/) { split($i, a, ":"); rq = a[3]; }
        }

        if (qs != "." && qe != "." && seq > 0) {
            diff = qe - qs - seq;
            score = 1 - (diff < 0 ? -diff : diff) / (qe - qs);
        } else {
            score = "NA";
        }

        split($1, id_parts, "/")
        current = id_parts[1]

        if (prev != "" && current != prev) {
            if (total > 0)
                printf "%-30s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", \
                    prev, total, lt60/total, lt50/total, lt40/total, lt30/total, rq_below/total
            total = lt60 = lt50 = lt40 = lt30 = rq_below = 0
        }

        prev = current
        total++
        if (mapq < 60) lt60++
        if (mapq < 50) lt50++
        if (mapq < 40) lt40++
        if (mapq < 30) lt30++
        if (rq != "." && rq < 0.99) rq_below++
    }
    END {
        if (total > 0)
            printf "%-30s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", \
                prev, total, lt60/total, lt50/total, lt40/total, lt30/total, rq_below/total
    }' >> "$output_count"

done < "$input_file"

