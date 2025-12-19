#!/bin/bash
ml 
in_dir=$1
kmer_size=$2
all_deduped_kmers=${in_dir}/deduped_kmer_counts.txt
not_included_kmers=${in_dir}/not_included_kmers.txt


# Match fasta files to variable names
while IFS=$'\t' read -r var_name var_value; do
    # Replace dots and hyphens with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    clean_var_name="${clean_var_name//-/_}"
    declare "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt

for file in ${in_dir}/*.ignored.bed; do
    # Generate tabular bed for ignored sequences [name, start, end, sequence]
    file_basename=$(basename "$file" .ignored.bed)
    clean_sample_id="${file_basename//\./_}"
    clean_sample_id="${clean_sample_id//-/_}"
    fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        echo "Processing ignored kmers for $file_basename"
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab >> $temp_bed
        rm "$temp_fasta"
    else
        echo "Processing ignored kmers for $file_basename"
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab >> $temp_bed
    fi
    #Extract kmers from bed file but keep only those not in deduped kmers + retain coordinates
    while IFS=$'\t' read -r seq_name start end sequence; do
        sequence_len=${#sequence}
        # Extract kmers of specified size from the sequence, search full kmer list.
        # If it's not there, add position of kmer to output file.
        for (( i=0; i<=sequence_len-kmer_size; i++ )); do
            kmer="${sequence:i:kmer_size+i-i}"
            kmer=$(echo "$kmer" | tr '[:lower:]' '[:upper:]')
            if ! grep -q -w "$kmer" "$all_deduped_kmers"; then
                kmer_start=$((start + i))
                kmer_end=$((start + i + kmer_size - 1))
                echo -e "${file_basename}\t${seq_name}\t${kmer_start}\t${kmer_end}\t${kmer}" >> ${not_included_kmers}
            fi
        done
    done < $temp_bed
done