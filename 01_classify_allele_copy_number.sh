#!/bin/bash

input_file="allele.adjusted.txt"

# Output files
out_4_no_na="4_no_na.txt"
out_1_na="1_na.txt"
out_2_na="2_na.txt"
out_3_na="3_na.txt"

# Clear the output files
> "$out_4_no_na"
> "$out_1_na"
> "$out_2_na"
> "$out_3_na"

awk -F'\t' '
NR==1 {
    # Write the header to all output files
    print > "'$out_4_no_na'"
    print > "'$out_1_na'"
    print > "'$out_2_na'"
    print > "'$out_3_na'"
    next
}
{
    # Check only columns 4–7 (Alleles A–D)
    na_count=0
    for(i=4;i<=7;i++){
        if($i=="NA") na_count++
    }
    if(na_count==0) print >> "'$out_4_no_na'"
    else if(na_count==1) print >> "'$out_1_na'"
    else if(na_count==2) print >> "'$out_2_na'"
    else if(na_count==3) print >> "'$out_3_na'"
}' "$input_file"

echo "  $out_4_no_na"
echo "  $out_1_na"
echo "  $out_2_na"
echo "  $out_3_na"
