# Loop through each unzipped BED file
for file in *.model.pbmm2.combined.bed; do
    bedtools intersect -a "$file" -b promoter.bed -wa -wb > "${file}.promoter.log"
done

output="qc_promoter_methylation_summary.tsv"
echo -e "File\tSum\tCount\tMedian\tBin1[0-20]\tBin2[20-40]\tBin3[40-60]\tBin4[60-80]\tBin5[80-100]" > "$output"

for file in *.model.pbmm2.combined.bed.promoter.log; do
  cut -f1-4 "$file" | sort -u | \
  gawk -v fname="$(basename "$file")" '
  {
    v = $4 + 0
    sum += v
    count++
    values[count] = v

    if (v >= 0 && v < 20) bin1++
    else if (v >= 20 && v < 40) bin2++
    else if (v >= 40 && v < 60) bin3++
    else if (v >= 60 && v < 80) bin4++
    else if (v >= 80 && v <= 100) bin5++
  }
  END {
    if (count == 0) {
      # output empty stats with zeros
      printf("%s\t0\t0\tNA\t0\t0\t0\t0\t0\n", fname)
      exit
    }
    asort(values)
    if (count % 2 == 1) {
      median = values[int((count+1)/2)]
    } else {
      median = (values[int(count/2)] + values[int(count/2)+1]) / 2
    }
    printf("%s\t%.0f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
           fname, sum, count, median,
           bin1/count, bin2/count, bin3/count, bin4/count, bin5/count)
  }'
done >> "$output"
