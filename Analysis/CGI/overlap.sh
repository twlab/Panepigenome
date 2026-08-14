#!/usr/bin/env bash
fea=CpG_shelf
type_tag=hap2
AWK=awk

while read -r sample_id; do
    input_txt="../${sample_id}.${type_tag}.ins.meth.bed"

    # Resolve matched hap naming
    if [[ "$type_tag" == "pat" ]]; then
        type_tag1="hap1"
    elif [[ "$type_tag" == "mat" ]]; then
        type_tag1="hap2"
    else
        type_tag1="$type_tag"
    fi

    matched="/scratch/zdong/Projects/PanEpiG/V1-9/PAV_variants/Phased/Merge/CGI/INS_"$fea"/CpG/${sample_id}.${type_tag1}.CpG.overlapped"

    out_overlap="${sample_id}.${type_tag}.del.$fea.meth.bed"

    echo "Processing: $sample_id $type_tag"

    # 清空输出文件
    : > "$out_overlap"

    # 按坐标匹配
    $AWK -F'\t' -v OFS="\t" -v out1="$out_overlap" '
    NR==FNR {
        # File B (matched): chr start end in columns 1,2,3; extra cols 10-15
        key = $1 FS $2 FS $3
        extra[key] = $10 FS $11 FS $12 FS $13 FS $14 FS $15
        next
    }
    {
        # File A (input): chr start end in columns 10,11,12
        key = $10 FS $11 FS $12
        if (key in extra)
            print $0, extra[key] > out1
    }
    ' "$matched" "$input_txt"

done < input.log
