#!/bin/bash

in_dir=$1
kmer_size=$2
mem=$3
threads=$4

output_deduped_file=${in_dir}/deduped_fasta_kmer_counts.txt
output_ignored_file=${in_dir}/ignored_fasta_kmer_counts.txt

echo ">kmer_counts" > ${in_dir}/all_sequences.txt
#Find all fasta files used in the dataset
sed 's/$/\n>/' ${in_dir}/all_dev.txt >> ${in_dir}/all_sequences.txt
sed '$!s/$/\n>/' ${in_dir}/all_train.txt >> ${in_dir}/all_sequences.txt


while IFS=$'\t' read -r var_name var_value; do
    # Replace dots with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    declare "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt


##Create sequnce file for ignored sequences
ignored_files=$(ls ${in_dir}/*.ignored.bed)
ignored_sequences=${in_dir}/all_ignored.txt
if [ -e $ignored_sequences ]; then
    rm $ignored_sequences
fi
touch $ignored_sequences

#For each file, getfasta, then pipe that txt to a bedfile with [chrom, start, end, sample_id, sequence]
for file in $ignored_files; do
    file_basename=$(basename "$file" .ignored.bed)
    clean_sample_id="${file_basename//\./_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -name >> $ignored_sequences
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -name >> $ignored_sequences
    fi
done

#Create a temporary file to store all sequences
mkdir -p ${in_dir}/kmc_dedup_temp
mkdir -p ${in_dir}/kmc_ignored_temp


kmc -k${kmer_size} -fm -ci1 -t${threads} -m${mem} ${in_dir}/all_sequences.txt ${in_dir}/deduped_fasta_files ${in_dir}/kmc_dedup_temp
kmc_tools transform ${in_dir}/deduped_fasta_files dump ${output_deduped_file}

kmc -k${kmer_size} -fm -ci1 -t${threads} -m${mem} ${ignored_sequences} ${in_dir}/ignored_fasta_files ${in_dir}/kmc_ignored_temp
kmc_tools transform ${in_dir}/ignored_fasta_files dump ${output_ignored_file}

total_ignored_and_deduped_file=${in_dir}/total_ignored_and_deduped_kmer_counts.txt
awk 'NR==FNR {sum[$1]+=$2; next} {sum[$1]+=$2} END {for (key in sum) print key, sum[key]}' ${output_deduped_file} ${output_ignored_file} > ${total_ignored_and_deduped_file}
rm -r ${in_dir}/*kmc_*