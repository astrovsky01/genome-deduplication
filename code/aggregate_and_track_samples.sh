#!/usr/bin/env bash

# Produces a corresponding txt file for each input bed
# e.g. for the above example, produces human.dev.txt and human.train.txt

# Detect OS and set appropriate decompression command
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    ZCAT_CMD="gzcat"
else
    # Linux and others
    ZCAT_CMD="zcat"
fi

# module load bedtools  # Comment out if not on HPC system

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

#For each file, getfasta, then pipe that txt to a bedfile with [chrom, start, end, sample_id, sequence]
for file in $train_files; do
    file_basename=$(basename "$file" .samples.train.bed)
    clean_sample_id="${file_basename//\./_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        $ZCAT_CMD "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab > $temp_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab > $temp_bed
    fi
    while IFS=$'\t' read -r coord sequence; do
        echo "$coord" | awk -F':|-' -v fasta_file="$fasta_file" -v sample_id="$file_basename" -v seq="$sequence" '{print $1 "\t" $2 "\t" $3 "\t" fasta_file "\t" sample_id "\t" seq}' >> $train_bed
    done < "$temp_bed"
    rm "$temp_bed"
done

for file in $dev_files; do
    file_basename=$(basename "$file" .samples.dev.bed)
    clean_sample_id="${file_basename//\./_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        $ZCAT_CMD "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab > $temp_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab > $temp_bed
    fi
    while IFS=$'\t' read -r coord sequence; do
        echo "$coord" | awk -F':|-' -v fasta_file="$fasta_file" -v sample_id="$file_basename" -v seq="$sequence" '{print $1 "\t" $2 "\t" $3 "\t" fasta_file "\t" sample_id "\t" seq}' >> $dev_bed
    done < "$temp_bed"
    rm "$temp_bed"
done

#Now shuffle the train and dev bed files
python ./code/shuffle.py $train_bed $random_seed ${train_bed}.shuf
python ./code/shuffle.py $dev_bed $random_seed ${dev_bed}.shuf
mv ${train_bed}.shuf $train_bed
mv ${dev_bed}.shuf $dev_bed

#Now separate columns 12,3,and 5 to bed files and 6 to txt file
cut -f6 $train_bed > ${train_txt}
cut -f6 $dev_bed > ${dev_txt}
cut -f1,2,3,5 $train_bed > ${train_bed}.tmp
mv ${train_bed}.tmp $train_bed
cut -f1,2,3,5 $dev_bed > ${dev_bed}.tmp
mv ${dev_bed}.tmp $dev_bed