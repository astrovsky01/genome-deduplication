#!/usr/bin/env bash

# Input arg
data_folder=tests/out/small_input

# Global args
dev_pct=0.1
max_dev_samples=100
shuf_seed=123

## Create train/dev txt files for each sample in the input directory

# Get all samples.bed files in the input directory
sample_beds=$(ls $data_folder/*.samples.bed)
for sample_bed in $sample_beds; do
	
	# Get all needed file names
	basename=$(basename $sample_bed | rev | cut -d'.' -f3- | rev)
	echo $basename
	sample_fasta=$(grep "^${basename}" ${data_folder}/basename_fasta_match.txt | cut -f2) # query corresponding fasta for this sample from basename_fasta_match.txt
	dev_sample_bed=${data_folder}/${basename}.samples.dev.bed
	train_sample_bed=${data_folder}/${basename}.samples.train.bed
	dev_sample_txt=${data_folder}/${basename}.samples.dev.txt
	train_sample_txt=${data_folder}/${basename}.samples.train.txt

	# Get training/dev beds/txts for each sample bed
	echo "Partitioning samples"
	./code/partition_samples.sh $sample_bed $dev_pct $max_dev_samples
	echo "Getting sequences"
	./code/get_sample_seqs.sh $sample_fasta $dev_sample_bed $train_sample_bed

done

## Group all sample train/dev files into global ones for the directory
#python code/aggregate_samples.py $data_folder $shuf_seed

