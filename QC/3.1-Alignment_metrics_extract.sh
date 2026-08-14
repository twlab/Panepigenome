#!/bin/bash

dir="/storage2/fs1/hprc/Active/zdong/V1-50/QC/step2-aligment_QC_106hg002hg005/Aligh_t2t"
output="summary_alignment_metrics.tsv"

# Write header
echo -e "Sample\tFile\tMapped_Reads\tAlignments\tMapped_Bases\tIdentity\tMax_Read_Length\tMean_Read_Length" > "$output"

# Process each alignment output log file in the dir
for log in "$dir"/*.log; do
  awk -v logfile="$log" '
  BEGIN {
      FS = ": ";
      OFS = "\t";
  }
  {
      if ($0 ~ /READ input file:/) {
          n = split($2, parts, "/");
          sample = parts[n-1];
          filename = parts[n];
      }
      if ($0 ~ /Mapped Reads:/) mapped_reads = $2;
      if ($0 ~ /Alignments:/) alignments = $2;
      if ($0 ~ /Mapped Bases:/) mapped_bases = $2;
      if ($0 ~ /Mean Gap-Compressed Sequence Identity:/) identity = $2;
      if ($0 ~ /Max Mapped Read Length:/) max_len = $2;
      if ($0 ~ /Mean Mapped Read Length:/) mean_len = $2;
  }
  END {
      if (sample != "" && filename != "")
          print sample, filename, mapped_reads, alignments, mapped_bases, identity, max_len, mean_len;
      else
          print "MISSING", logfile, "-", "-", "-", "-", "-", "-";
  }' "$log" >> "$output"
done

