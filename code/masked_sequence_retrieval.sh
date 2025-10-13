#!/usr/bin/env bash



mask=false
ignore=false

while getopts ":i:mn" opt; do
    case $opt in
        i)
            data_folder="$OPTARG"
            ;;
        m)
            mask=true
            ;;
        n)
            ignore=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i input_folder [-m] [-n]"
            echo "  -m : process .masks.bed files to get masked sequences"
            echo "  -n : process .ignored.bed files to get ignored sequences"
            exit 1
            ;;

    esac
done

in_dir=$data_folder

if $mask; then
    masked_beds=$(ls $data_folder/*.masks.bed)
fi

if $ignore; then
    ignored_beds=$(ls $data_folder/*.ignored.bed)
fi

map_file=${in_dir}/basename_fasta_match.txt
while IFS=$'\t' read -r var_name var_value; do
    # Replace dots with underscores to make valid bash variable names
    clean_var_name="${var_name//\./_}"
    declare "$clean_var_name=$var_value"
    echo "$clean_var_name=$var_value"
done < ${in_dir}/basename_fasta_match.txt

if $mask; then
    for bed in $masked_beds; do
        file_basename=$(basename $bed .masks.bed)
        clean_sample_id="${file_basename//\./_}"
        fasta_file="${!clean_sample_id}"  # Get actual fasta file from nam
        outfile=${bed%.masks.bed}.masks.fa
        if [[ "$fasta_file" == *.gz ]]; then
            temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
            gunzip -c "$fasta_file" > "$temp_fasta"
            bedtools getfasta -fi "$temp_fasta" -bed "$bed" -fo  $outfile
            rm "$temp_fasta"
        else
            bedtools getfasta -fi "$fasta_file" -bed "$bed" -fo $outfile
        fi
    done
fi

if $ignore; then
    for bed in $ignored_beds; do
        file_basename=$(basename $bed .ignored.bed)
        clean_sample_id="${file_basename//\./_}"
        fasta_file="${!clean_sample_id}"  # Get actual fasta file from name
        outfile=${bed%.ignored.bed}.ignored.fa
        if [[ "$fasta_file" == *.gz ]]; then
            temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
            if [[ "$OSTYPE" == "darwin"* ]]; then
                gzcat "$fasta_file" > "$temp_fasta"
            else
                zcat "$fasta_file" > "$temp_fasta"
            fi
            bedtools getfasta -fi "$temp_fasta" -bed "$bed" -fo  $outfile
            rm "$temp_fasta"
        else
            bedtools getfasta -fi "$fasta_file" -bed "$bed" -fo $outfile
        fi
    done
fi