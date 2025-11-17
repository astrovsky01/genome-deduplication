#!/bin/bash

input_file=$1
base_name=$(basename "$input_file" .bed)
output_file=${base_name}_chunk_lengths.csv
awk '{printf "%s,", ($3 - $2)}' "$input_file" >> "$output_file"
sed -i 's/,$//' "$output_file"