#!/bin/bash
ml bedtools

in_dir=$1
all_deduped_kmers=$2
kmer_size=$3
threads=$4
mem=$5
not_included_kmers=${in_dir}/not_included_kmers.txt
script_dir=$(dirname "$0")
unused_sequence_txt=${in_dir}/unused_sequences.txt
unused_sequences_fasta=${in_dir}/unused_sequences.fasta
unused_kmers_txt=${in_dir}/unused_kmers.txt
input_file=${in_dir}/python_inputs.txt

if [ -e "$not_included_kmers" ]; then
    rm "$not_included_kmers"
fi
touch $not_included_kmers

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
    echo -e "$temp_bed\t$file_basename" >> $input_file
done


# python ${script_dir}/unused_kmers.py --ignored_file $temp_bed --file_name $file_basename --kmer_size $kmer_size --deduped_kmers $all_deduped_kmers --outfile $not_included_kmers

python ${script_dir}/unused_kmers.py --input $input_file --kmer_size $kmer_size --deduped_kmers $all_deduped_kmers --outfile $not_included_kmers

tmp_beds=$(cut -f1 $input_file)
rm $input_file
for bed in $tmp_beds; do
    rm $bed
done

# Break sequences into kmers and kmc to count unique unused kmers
mkdir -p ${in_dir}/kmc_tmp
cut -f 4 $not_included_kmers > $unused_sequence_txt
cut -f 1,2,3,5 $not_included_kmers > ${not_included_kmers%.txt}.bed

# Convert raw sequences to FASTA format for KMC
awk '{print ">"NR"\n"$0}' $unused_sequence_txt > $unused_sequences_fasta
kmc -k${kmer_size} -cs4294967295 -fm -b -ci1 -t${threads} -m${mem} $unused_sequences_fasta ${in_dir}/unused_kmers_db ${in_dir}/kmc_tmp
kmc_tools transform ${in_dir}/unused_kmers_db dump $unused_kmers_txt
rm ${in_dir}/unused_kmers_db*
# rm -rf ${in_dir}/kmc_tmp
# rm $unused_sequence_txt
# rm $unused_sequences_fasta
# rm $not_included_kmers