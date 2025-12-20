#!/bin/bash
ml bedtools

in_dir=$1
all_deduped_kmers=$2
kmer_size=$3
not_included_kmers=${in_dir}/not_included_kmers.txt
script_dir=$(dirname "$0")

if [ -e "$not_included_kmers" ]; then
    rm "$not_included_kmers"
fi
touch $not_included_kmers

# Match fasta files to variable names
while IFS=$'\t' read -r var_name var_value; do
    # Replace dots and hyphens with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    clean_var_name="${clean_var_name//-/_}"
    declare "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt

for file in ${in_dir}/*.ignored.bed; do
    # Generate tabular bed for ignored sequences [name, start, end, sequence]
    file_basename=$(basename "$file" .ignored.bed)
    clean_sample_id="${file_basename//\./_}"
    clean_sample_id="${clean_sample_id//-/_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        echo "Processing ignored kmers for $file_basename"
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab >> $temp_bed
        rm "$temp_fasta"
    else
        echo "Processing ignored kmers for $file_basename"
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab >> $temp_bed
    fi
    python ${script_dir}/unused_kmers.py --ignored_file $temp_bed --file_name $file_basename --kmer_size $kmer_size --deduped_kmers $all_deduped_kmers --outfile $not_included_kmers
    rm $temp_bed
done