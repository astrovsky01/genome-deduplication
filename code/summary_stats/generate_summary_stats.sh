#!/bin/bash

#Generate kmer counts

input_dir=$1
output_dir=${input_dir}/$2
kmer_size=$3
mem=$4
threads=$5

mkdir -p $output_dir
if [ -e ${output_dir}/summary_stats.txt ]; then
    rm ${output_dir}/summary_stats.txt
fi
touch ${output_dir}/summary_stats.txt

source_kmer_counts="/source_fasta_kmer_counts.txt"
deduped_file="/deduped_kmer_counts.txt"
ignored_file="/ignored_kmer_counts.txt"
ignored_and_deduped_file="/ignored_deduped_kmer_counts.txt"
masked_file="/masked_kmer_counts.txt"
combined_file="/combined_all_kmers.txt"


#generate dedup kmer counts
sh code/summary_stats/kmer_counter_dedup_files.sh $input_dir $kmer_size $mem $threads
mv ${input_dir}${deduped_file} ${output_dir}${deduped_file}
mv ${input_dir}${ignored_file} ${output_dir}${ignored_file}
mv ${input_dir}${masked_file} ${output_dir}${masked_file}
mv ${input_dir}${combined_file} ${output_dir}${combined_file}
mv ${input_dir}${ignored_and_deduped_file} ${output_dir}${ignored_and_deduped_file}

#Generate source data kmer counts
sh code/summary_stats/kmer_counter_original_files.sh $input_dir $kmer_size $mem $threads
mv ${input_dir}${source_kmer_counts} ${output_dir}${source_kmer_counts}

mv ${input_dir}/*.log ${output_dir}/

#Generate kmer count statistics
echo "Calculating summary statistics on source files"
unique_source_kmers=$(grep "No. of unique k-mers" ${output_dir}/source_kmer_counts.log | awk '{print $6}')
source_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/source_kmer_counts.log | awk '{print $6}')

echo "Calculating summary statistics on dedup files"
unique_deduped_kmers=$(grep "No. of unique k-mers" ${output_dir}/deduped_kmer_counts.log | awk '{print $6}')
deduped_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/deduped_kmer_counts.log | awk '{print $6}')

echo "Calculating summary statistics on ignored files"
unique_ignored_kmers=$(grep "No. of unique k-mers" ${output_dir}/ignored_kmer_counts.log | awk '{print $6}')
ignored_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/ignored_kmer_counts.log | awk '{print $6}')

echo "Calculating summary statistics on ignored + deduped files"
unique_ignored_deduped_kmers=$(grep "No. of unique k-mers" ${output_dir}/ignored_deduped_kmer_counts.log | awk '{print $6}')
deduped_ignored_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/ignored_deduped_kmer_counts.log | awk '{print $6}')

echo "Calculating summary statistics on masked files"
unique_masked_kmers=$(grep "No. of unique k-mers" ${output_dir}/masked_kmer_counts.log | awk '{print $6}')
masked_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/masked_kmer_counts.log | awk '{print $6}')

echo "Calculating summary statistics on combined files"
unique_combined_kmers=$(grep "No. of unique k-mers" ${output_dir}/combined_kmer_counts.log | awk '{print $6}')
combined_total_kmer_count=$(grep "Total no. of k-mers" ${output_dir}/combined_kmer_counts.log | awk '{print $6}')

# N Counter
echo "Counting ambiguous bases in source fasta files"
total=0
validchars="ACGTacgt"
while IFS= read -r file; do
    if [[ "$file" == *.gz ]]; then
        count=$(zgrep -v "^>" "$file" | grep -o "[^$validchars]" | wc -l)
    else
        count=$(grep -v "^>" "$file" | grep -o "[^$validchars]" | wc -l)
    fi
    total=$((total + count))
done < ${input_dir}/fasta_file_list.txt

echo "Counting ambigous bases found in dedpuplication"
ambiguous=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ${input_dir}/*.ambiguous.bed)

#Write summary stats to file
echo -e "Summary Statistics for kmer size ${kmer_size}:\n" >> ${output_dir}/summary_stats.txt
echo -e "Source\t|\tUnique kmer counts:\t|\tTotal kmer counts:\n" >> ${output_dir}/summary_stats.txt
echo -e "-----------------------------------------------\n" >> ${output_dir}/summary_stats.txt
echo -e "Original fasta files\t|\t${unique_source_kmers}\t|\t${source_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Deduplicated files\t|\t${unique_deduped_kmers}\t|\t${deduped_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Ignored files\t|\t${unique_ignored_kmers}\t|\t${ignored_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Ignored and deduplicated files\t|\t${unique_ignored_deduped_kmers}\t|\t${deduped_ignored_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Masked files\t|\t${unique_masked_kmers}\t|\t${masked_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Combined files\t|\t${unique_combined_kmers}\t|\t${combined_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "Total N bases in source fasta files: ${total}\n" >> ${output_dir}/summary_stats.txt
echo -e "Total Ambiguous bases found: ${ambiguous}\n" >> ${output_dir}/summary_stats.txt

echo "Generating spacing figures"
if [ ! -e ${output_dir}/spacing_figures ]; then
    mkdir ${output_dir}/spacing_figures
fi
while IFS=$'\t' read -r basename fasta_loc; do
    # Replace dots with underscores to make valid bash variable names
    sh code/summary_stats/calculate_distance_between_dedups.sh ${input_dir}/${basename}.fai ${input_dir}/${basename}.samples.bed ${input_dir}/${basename}_distance_between_deduplicated_regions.csv
    python code/summary_stats/spacing_info.py --input_file ${input_dir}/${basename}_distance_between_deduplicated_regions.csv --output_dir ${output_dir}/spacing_figures --type dedup --bincount 1000
    # rm ${input_dir}/${basename}_distance_between_deduplicated_regions.csv
    sh code/summary_stats/chunk_lengths.sh ${input_dir}/${basename}.masks.bed
    python code/summary_stats/spacing_info.py --input_file ${input_dir}/${basename}.masks_chunk_lengths.csv --output_dir ${output_dir}/spacing_figures --type masked --bincount 1000
    # rm ${input_dir}/${basename}.masks_chunk_lengths.csv
    sh code/summary_stats/chunk_lengths.sh ${input_dir}/${basename}.ignored.bed
    python code/summary_stats/spacing_info.py --input_file ${input_dir}/${basename}.ignored_chunk_lengths.csv --output_dir ${output_dir}/spacing_figures --type ignored --bincount 1000
    # rm ${input_dir}/${basename}.ignored_chunk_lengths.csv
done < ${input_dir}/basename_fasta_match.txt