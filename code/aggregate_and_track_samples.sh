#!/usr/bin/env bash

# Produces a corresponding txt file for each input bed
# e.g. for the above example, produces human.dev.txt and human.train.txt


# module load bedtools  # Comment out if not on HPC system

in_dir=$1
random_seed=$2
test_pct=${3:-0}
no_txt=${4:-false}
train_bed=${in_dir}/all_train.bed
dev_bed=${in_dir}/all_dev.bed
train_txt=${in_dir}/all_train.txt
dev_txt=${in_dir}/all_dev.txt
touch $train_bed
touch $dev_bed
if [[ $no_txt == false ]]; then
    touch $train_txt
    touch $dev_txt
fi
if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
    test_bed=${in_dir}/all_test.bed
    test_txt=${in_dir}/all_test.txt
    touch $test_bed
    test_files=$(ls ${in_dir}/*.samples.test.bed)
    if [[ $no_txt == false ]]; then
        touch $test_txt
    fi
fi
train_files=$(ls ${in_dir}/*.samples.train.bed)
dev_files=$(ls ${in_dir}/*.samples.dev.bed)

#For each file, getfasta, then pipe that txt to a bedfile with [chrom, start, end, sample_id, sequence]
for file in $train_files; do
    file_basename=$(basename "$file" .samples.train.bed)
    # For shuffle, also remove the .final suffix if present
    file_basename=$(basename "$file_basename" .final)
    fasta_file=$(grep $file_basename ${in_dir}/basename_fasta_match.txt | cut -f2)
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab > $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $train_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $train_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab > $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $train_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $train_bed
    fi
    rm $temp_bed
done

for file in $dev_files; do
    file_basename=$(basename "$file" .samples.dev.bed)
    file_basename=$(basename "$file_basename" .final)
    fasta_file=$(grep $file_basename ${in_dir}/basename_fasta_match.txt | cut -f2)
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab >> $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $dev_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $dev_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab >> $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $dev_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $dev_bed
    fi
    rm $temp_bed
done

if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
    for file in $test_files; do
        file_basename=$(basename "$file" .samples.test.bed)
        file_basename=$(basename "$file_basename" .final)
        fasta_file=$(grep $file_basename ${in_dir}/basename_fasta_match.txt | cut -f2)
        temp_bed=$(mktemp)
        if [[ "$fasta_file" == *.gz ]]; then
            temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
            gunzip -c "$fasta_file" > "$temp_fasta"
            bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab >> $temp_bed
            if [ "$no_txt" = false ]; then
                cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $test_txt
            fi
            awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $test_bed
            rm "$temp_fasta"
        else
            bedtools getfasta -fi "$fasta_file" -bed "$file" -tab >> $temp_bed
            if [ "$no_txt" = false ]; then
                cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $test_txt
            fi
            awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $test_bed
        fi
        rm $temp_bed
    done
fi

#Now shuffle the train and dev bed files
shuf --random-source="$train_bed" -o ${train_bed}.shuf $train_bed
if [ "$no_txt" = false ]; then
    shuf --random-source="$train_bed" -o ${train_txt}.shuf $train_txt
fi
shuf --random-source="$dev_bed" -o ${dev_bed}.shuf $dev_bed
if [ "$no_txt" = false ]; then
    shuf --random-source="$dev_bed" -o ${dev_txt}.shuf $dev_txt
fi
if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
    shuf --random-source="$test_bed" -o ${test_bed}.shuf $test_bed
    if [ "$no_txt" = false ]; then
        shuf --random-source="$test_bed" -o ${test_txt}.shuf $test_txt
    fi
fi
mv ${train_bed}.shuf $train_bed
mv ${dev_bed}.shuf $dev_bed
if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
    mv ${test_bed}.shuf $test_bed
fi
if [ "$no_txt" = false ]; then
    mv ${train_txt}.shuf $train_txt
    mv ${dev_txt}.shuf $dev_txt
    if awk -v pct="$test_pct" 'BEGIN { exit !(pct > 0) }'; then
        mv ${test_txt}.shuf $test_txt
    fi
fi
