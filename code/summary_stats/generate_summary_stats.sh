#!/bin/bash

#Generate kmer counts

input_dir=$1
output_dir=${input_dir}/$2
kmer_size=$3
mem=$4
threads=$5

mkdir -p $output_dir
touch ${output_dir}/summary_stats.txt

source_kmer_counts="/source_data_kmer_counts.txt"
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
sh code/summary_stats/kmer_counter_source_data.sh $input_dir $kmer_size $mem $threads
mv ${input_dir}${source_kmer_counts} ${output_dir}${source_kmer_counts}


#Generate kmer count statistics
unique_source_kmers=$(wc -l ${output_dir}${source_kmer_counts})
source_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${source_kmer_counts})

unique_deduped_kmers=$(wc -l ${output_dir}${deduped_file})
deduped_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${deduped_file})

unique_ignored_kmers=$(wc -l ${output_dir}${ignored_file})
ignored_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${ignored_file})

unique_ignored_deduped_kmers=$(wc -l ${output_dir}${ignored_and_deduped_file})
deduped_ignored_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${ignored_and_deduped_file})

unique_masked_kmers=$(wc -l ${output_dir}${masked_file})
masked_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${masked_file})

unique_combined_kmers=$(wc -l ${output_dir}${combined_file})
combined_total_kmer_count=$(awk '{sum+=$2} END{print sum}' ${output_dir}${combined_file})


#Write summary stats to file
echo -e "Summary Statistics for kmer size ${kmer_size}:\n" >> ${output_dir}/summary_stats.txt
echo -e "Source\tUnique kmer counts:\tTotal kmer counts:\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_source_kmers}\t${source_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_deduped_kmers}\t${deduped_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_ignored_kmers}\t${ignored_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_ignored_deduped_kmers}\t${deduped_ignored_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_masked_kmers}\t${masked_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt
echo -e "${unique_combined_kmers}\t${combined_total_kmer_count}\n" >> ${output_dir}/summary_stats.txt



#Show kmer counts of original dataset overlayed with deduped and deduped+ignored kmer distributions
python code/summary_stats/kmer_count_summary.py ${output_dir}

#Show where kmers came from in deduped datasets
python code/summary_stats/kmer_distribution_plot.py ${input_dir} ${output_dir}