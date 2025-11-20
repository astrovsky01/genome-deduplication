#!/bin/bash

input_file=$1
base_name=$(basename "$input_file" .bed)
output_dir=$(dirname "$input_file")
output_file=${output_dir}/${base_name}_chunk_lengths.csv
if [ -e "$output_file" ]; then
    rm "$output_file"
fi
touch "$output_file"
awk '{printf "%s,", ($3 - $2)}' "$input_file" >> "$output_file"
sed -i 's/,$//' "$output_file"