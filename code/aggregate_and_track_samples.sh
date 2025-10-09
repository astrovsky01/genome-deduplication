#!/usr/bin/env bash

# Usage: ./code/get_sample_seqs.sh [fasta] [beds...]
# e.g.:  ./code/get_sample_seqs.sh human.fa human.dev.bed human.train.bed
# Produces a corresponding txt file for each input bed
# e.g. for the above example, produces human.dev.txt and human.train.txt

module load bedtools

#in_bed=tests/out/dummy/dummy_1.samples.bed
#in_fa=tests/test-data/dummy_1.fa
#out_file=foo.txt

#in_bed=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.bed
#in_fa=/scratch4/mschatz1/mrevsin1/ncbi_genomes/genomes/GCA_001775245.1_ASM177524v1_genomic.fna.gz
#out_file=tests/out_small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.txt

in_dir=$1
random_seed=$2

train_bed=${in_dir}/all_train.bed
dev_bed=${in_dir}/all_dev.bed
train_txt=${in_dir}/all_train.txt
dev_txt=${in_dir}/all_dev.txt
touch $train_bed
touch $dev_bed

train_files=$(ls ${in_dir}/*.samples.train.bed)
dev_files=$(ls ${in_dir}/*.samples.dev.bed)

while IFS=$'\t' read -r var_name var_value; do
    # Replace dots with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    declare "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt


# Create combined train and test bed files with labeled origins
for file in $train_files; do
    file_basename=$(basename "$file" .samples.train.bed)
    cat "$file" | awk -v fname="$file_basename" '{print $0 "\t" fname}' >> $train_bed
done
python code/shuffle.py $train_bed $random_seed $train_bed.shuf
for file in $dev_files; do
    file_basename=$(basename "$file" .samples.dev.bed)
    cat "$file" | awk -v fname="$file_basename" '{print $0 "\t" fname}' >> $dev_bed
done
python code/shuffle.py $dev_bed $random_seed $dev_bed.shuf

while IFS=$'\t' read -r chr start end sample_id; do
    # Convert dots to underscores to match declared variable names
    clean_sample_id="${sample_id//\./_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    if [[ "$fasta_file" == *.gz ]]; then
        echo -e "$chr\t$start\t$end" | bedtools getfasta -fi <(zcat "$fasta_file") -bed - -tab | cut -f2
    else
        echo -e "$chr\t$start\t$end" | bedtools getfasta -fi "$fasta_file" -bed - -tab | cut -f2
    fi
done < "$train_bed" >> "${train_txt}"


while IFS=$'\t' read -r chr start end sample_id; do
    # Convert dots to underscores to match declared variable names
    clean_sample_id="${sample_id//\./_}"
    fasta_file="${!clean_sample_id}"
    if [[ "$fasta_file" == *.gz ]]; then
        echo -e "$chr\t$start\t$end" | bedtools getfasta -fi <(zcat "$fasta_file") -bed - -tab | cut -f2
    else
        echo -e "$chr\t$start\t$end" | bedtools getfasta -fi "$fasta_file" -bed - -tab | cut -f2
    fi
done < "$dev_bed" >> "${dev_txt}"
