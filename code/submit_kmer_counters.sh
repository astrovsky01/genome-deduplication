#!/bin/bash

#SET SUBMISSION VALUES FOR BOTH SCRIPTS
# =======================
in_dir=$1
kmer_size=$2
mem=$3
threads=$4
job_original="base_fasta_files"
job_dedup="deduped_fasta_files"
allocation="mschatz1"
partition="parallel"
time="1:00:00"
root_dir="/scratch4/mschatz1/aostrov4/retraining_project/genome-deduplication"
# =======================

echo "Submitting ORIGINAL kmer counting job..."
#SUBMIT KMER COUNTER FOR ORIGINAL FILES
sbatch \
    --mem=${mem}G \
    --ntasks-per-node=${threads} \
    --job-name=kmer_original \
    -A ${allocation} \
    --output=${in_dir}/${job_original}.out  \
    --error=${in_dir}/${job_original}.err \
    --partition=${partition} \
    --time=${time} \
    --wrap="cd ${root_dir} && sh code/kmer_counter_original_files.sh ${in_dir} ${kmer_size} ${mem} ${threads}"

echo "Submitting DEDUPED kmer counting job..."
#SUBMIT KMER COUNTER FOR DEDUPED FILES
sbatch \
    --mem=${mem}G \
    --ntasks-per-node=${threads} \
    --job-name=kmer_dedup \
    -A ${allocation} \
    --output=${in_dir}/${job_dedup}.out  \
    --error=${in_dir}/${job_dedup}.err \
    --partition=${partition} \
    --time=${time} \
    --wrap="cd ${root_dir} && sh code/kmer_counter_dedup_files.sh ${in_dir} ${kmer_size} ${mem} ${threads}"