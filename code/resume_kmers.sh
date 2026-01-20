#!/bin/bash
ml bedtools

in_dir=$1
all_deduped_kmers=$2
kmer_size=$3
threads=$4
mem=$5
previous_kmers=${in_dir}/previous_kmers.txt
script_dir=$(dirname "$0")
previous_sequence_txt=${in_dir}/previous_sequences.txt
previous_sequences_fasta=${in_dir}/previous_sequences.fasta
previous_kmers_txt=${in_dir}/previous_kmers.txt
input_file=${in_dir}/python_inputs.txt

if [ -e "$previous_kmers" ]; then
    rm "$previous_kmers"
fi
touch $previous_kmers

# Match fasta files to variable names
while IFS=$'\t' read -r var_name var_value; do
    # Replace dots and hyphens with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    clean_var_name="${clean_var_name//-/_}"
    declare "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt

for file in ${in_dir}/*.samples.bed; do
    # Generate tabular bed for ignored sequences [name, start, end, sequence]
    file_basename=$(basename "$file" .samples.bed)
    clean_sample_id="${file_basename//\./_}"
    clean_sample_id="${clean_sample_id//-/_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        echo "Processing sample kmers for $file_basename"
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" >> $sample_fasta
        rm "$temp_fasta"
    else
        echo "Processing sample kmers for $file_basename"
        bedtools getfasta -fi "$fasta_file" -bed "$file" >> $sample_fasta
    fi
done

kmc -k${kmer_size} -cs4294967295 -fm -b -ci1 -t${threads} -m${mem} $sample_fasta ${in_dir}/previous_kmers_db ${in_dir}/kmc_tmp
kmc_tools transform ${in_dir}/previous_kmers_db dump $previous_kmers_txt
rm ${in_dir}/previous_kmers_db*
rm $sample_fasta