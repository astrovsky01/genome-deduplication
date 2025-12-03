#!/bin/bash

#Go line-by-line for each contig
fai_file=$1
sample_file=$2
output_file=$3
if [ -e $output_file ]; then
    rm $output_file
fi
touch $output_file
if [ -e sample_contig.txt ]; then
    rm sample_contig.txt
fi
touch sample_contig.txt

while IFS=$'\t' read -r contig length offset linebases linewidth qualoffset
do
    grep "^${contig}	" $sample_file > sample_contig.txt
    first_start=$(head -1 sample_contig.txt | cut -f2)
    if [[ -n "$first_start" ]]; then
        # Output gap from start of contig (position 0) to first region
        printf "%s," "$first_start" >> $output_file
    fi
    
    distances=$(awk 'NR > 1 { 
        gap = $2 - prev_col3; 
        results = results "," gap 
    } { 
        prev_col3 = $3 
    } END { 
        if (results != "") printf "%s", substr(results, 2) 
    }' sample_contig.txt)
    if [[ -n "$distances" ]]; then
        printf "%s," "$distances" >> $output_file
    fi
    
    last_end=$(tail -1 sample_contig.txt | cut -f3)
    if [[ -n "$last_end" ]]; then
        remaining_length=$((length - last_end))
        printf "%s," "${remaining_length}" >> $output_file
    fi
done < $fai_file

# Remove the final comma
sed -i 's/,$//' $output_file
rm sample_contig.txt