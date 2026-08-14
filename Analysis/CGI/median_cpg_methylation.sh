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

