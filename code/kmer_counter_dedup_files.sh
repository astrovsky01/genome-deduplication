#!/bin/bash

in_dir=$1
kmer_size=$2
mem=$3
threads=$4

output_file=${in_dir}/deduped_fasta_kmer_counts.txt

echo ">kmer_counts" > ${in_dir}/all_sequences.txt
#Find all fasta files used in the dataset
sed 's/$/\n>/' ${in_dir}/all_dev.txt >> ${in_dir}/all_sequences.txt
sed '$!s/$/\n>/' ${in_dir}/all_train.txt >> ${in_dir}/all_sequences.txt


#Create a temporary file to store all sequences
mkdir -p ${in_dir}/kmc_dedup_temp


kmc -k${kmer_size} -fm -ci1 -t${threads} -m${mem} ${in_dir}/all_sequences.txt deduped_fasta_files ${in_dir}/kmc_dedup_temp
kmc_tools transform deduped_fasta_files dump ${output_file}
rm -r ${in_dir}/kmc_dedup_temp
rm ${in_dir}/all_sequences.txt
echo "" >> ${in_dir}/all_sequences.txt