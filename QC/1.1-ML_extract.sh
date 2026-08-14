while read -r bam; do
  samtools view "$bam" -@8 | awk -v bamname="$bam" '
    {
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^[Mm][Ll]:B:C,/) {
          ml = $i
          sub(/^[Mm][Ll]:B:C,/, "", ml)
          n = split(ml, vals, ",")
          for (j = 1; j <= n; j++) {
            val = vals[j]
            if (val >= 1 && val <= 255) {
              # Divide into 5 bins, each 51 units
              bin = int((val - 1) / 51) + 1
              bins[bin]++
              sum += val
              count++
            }
          }
        }
      }
    }
    END {
      if (count > 0) {
        printf "%s\t%d\t%d\t%.4f", bamname, sum, count, sum / count
        for (b = 1; b <= 5; b++) {
          ratio = (bins[b] > 0) ? bins[b] / count : 0
          printf "\t%.4f", ratio
        }
        printf "\n"
      } else {
        printf "%s\tNo ML values found.\n", bamname
      }
    }'
done < example_input_1.1-ML_extract.log
