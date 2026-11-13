#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Usage: ./code/partition_samples.sh [X.samples.bed] [dev set %] [max dev set count] [random seed] [test set %]
# e.g.:  ./code/partition_samples.sh human.samples.bed 0.1 100000 123 0.1
# Produces output dev and train bed files. If test set percentage is greater than zero, also produces a test bed file.

# Input args, hard-coded for now
#sample_bed=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.bed
#dev_pct=0.1 # Really dev fraction
#max_dev_samples=100
sample_bed=$1
dev_pct=$2
max_dev_samples=$3
random_seed=$4
test_pct=$5


# Get number of samples in this file
n_samples=$(wc -l $sample_bed)
n_samples=$(echo $n_samples | cut -d' ' -f1)
# Get the number of dev samples based on the dev_pct
n_dev_samples=$(echo "$n_samples * $dev_pct" | bc -l)
n_dev_samples=$(echo "$n_dev_samples/1" | bc)
n_test_samples=$(echo "$n_samples * $test_pct" | bc -l)
n_test_samples=$(echo "$n_test_samples/1" | bc)
# Make sure the dev sample count doesn't exceed the max dev samples count
if [[ $n_dev_samples -gt $max_dev_samples ]]; then
	n_dev_samples=$max_dev_samples
fi

# Get filename up to the suffix (.samples.bed)
file_basename="${sample_bed%.samples.bed}"
# Shuffle file to temporary file
shuf_file=${file_basename}.shuf.samples.bed
python "$SCRIPT_DIR/shuffle.py" "$sample_bed" "$random_seed" "$shuf_file"
#shuf $sample_bed > $shuf_file

# Partition shuffled sample file
if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
    head -n $n_test_samples $shuf_file > ${file_basename}.samples.test.bed
    head -n $((n_test_samples+n_dev_samples)) $shuf_file \
        | tail -n $n_dev_samples > ${file_basename}.samples.dev.bed
    tail -n +$((n_test_samples+n_dev_samples+1)) $shuf_file > ${file_basename}.samples.train.bed
else
    head -n $n_dev_samples $shuf_file > ${file_basename}.samples.dev.bed
    tail -n +$((n_dev_samples+1)) $shuf_file > ${file_basename}.samples.train.bed
fi

# Remove temporary file
rm -f $shuf_file
