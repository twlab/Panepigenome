file=matched_meth.bed.filtered
awk '{if(!/^#/ && $4!="") print $4}' $file | sort -n | awk '
{
  a[NR]=$1
}
END {
  if (NR % 2) 
    print a[(NR + 1) / 2]
  else
    print (a[NR / 2] + a[(NR / 2) + 1]) / 2
}'
