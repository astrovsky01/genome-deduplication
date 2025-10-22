#!/bin/bash

in_dir=$1
kmer_size=$2
mem=$3
threads=$4

output_file=${in_dir}/source_fasta_kmer_counts.txt
touch ${in_dir}/fasta_file_list.txt
#Find all fasta files used in the dataset
cut -f2 ${in_dir}/basename_fasta_match.txt > ${in_dir}/fasta_file_list.txt


#Create a temporary file to store all sequences

mkdir -p ${in_dir}/kmc_temp

kmc -k${kmer_size} -b -fm -ci1 -t${threads} -m${mem} @${in_dir}/fasta_file_list.txt original_fasta_files ${in_dir}/kmc_temp
kmc_tools transform original_fasta_files dump ${output_file}
rm -r ${in_dir}/kmc_temp