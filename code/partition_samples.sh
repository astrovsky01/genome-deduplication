#!/usr/bin/env bash

dev_pct=0.1
max_dev_samples=100
sample_file=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.bed

# Get number of samples in this file
n_samples=$(wc -l $sample_file | cut -d' ' -f1)

# Get the number of dev samples based on the dev_pct
n_dev_samples=$(echo "$n_samples * $dev_pct" | bc -l)
n_dev_samples=$(echo "$n_dev_samples/1" | bc)

# Make sure the dev sample count doesn't exceed the max dev samples count
if [[ $n_dev_samples -gt $max_dev_samples ]]; then
	n_dev_samples=$max_dev_samples
fi

# Get filename up to the suffix (.samples.bed)
file_basename=$(echo "${sample_file::-12}")

# Shuffle file to temporary file
shuf_file=${file_basename}.shuf.samples.bed
shuf $sample_file > $shuf_file

# Partition shuffled sample file
head -n $n_dev_samples $shuf_file > ${file_basename}.samples.dev.bed
tail -n +$((n_dev_samples+1)) $shuf_file > ${file_basename}.samples.train.bed

# Remove temporary file
rm -f $shuf_file

