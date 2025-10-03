#!/usr/bin/env bash

# Usage: ./code/partition_samples.sh [X.samples.bed] [dev set %] [max dev set count] [random seed]
# e.g.:  ./code/partition_samples.sh human.samples.bed 0.1 100000 123
# Produces an output dev and train bed file at X.samples.dev.bed and X.samples.train.bed

# Input args, hard-coded for now
#sample_bed=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.bed
#dev_pct=0.1 # Really dev fraction
#max_dev_samples=100
sample_bed=$1
dev_pct=$2
max_dev_samples=$3
random_seed=$4

# Get number of samples in this file
n_samples=$(wc -l $sample_bed | cut -d' ' -f1)

# Get the number of dev samples based on the dev_pct
n_dev_samples=$(echo "$n_samples * $dev_pct" | bc -l)
n_dev_samples=$(echo "$n_dev_samples/1" | bc)

# Make sure the dev sample count doesn't exceed the max dev samples count
if [[ $n_dev_samples -gt $max_dev_samples ]]; then
	n_dev_samples=$max_dev_samples
fi

# Get filename up to the suffix (.samples.bed)
file_basename=$(echo "${sample_bed::-12}")

# Shuffle file to temporary file
shuf_file=${file_basename}.shuf.samples.bed
python code/shuffle.py $sample_bed $random_seed $shuf_file
#shuf $sample_bed > $shuf_file

# Partition shuffled sample file
head -n $n_dev_samples $shuf_file > ${file_basename}.samples.dev.bed
tail -n +$((n_dev_samples+1)) $shuf_file > ${file_basename}.samples.train.bed

# Remove temporary file
rm -f $shuf_file

