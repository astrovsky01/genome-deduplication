#!/usr/bin/env bash


# Input arg
data_folder=$1

# Global args
dev_pct=0.1
max_dev_samples=100
shuf_seed=123


while getopts ":i:p:m:s:" opt; do
    case $opt in
        i)
            data_folder="$OPTARG"
            ;;
        p)
            dev_pct="$OPTARG"
            ;;
        m)
            max_dev_samples="$OPTARG"
            ;;
        s)
            shuf_seed="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i input_folder -p output_folder -m max_dev_samples -s seed"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo "Usage: $0 -i input_folder -p output_folder -m max_dev_samples -s seed"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided and not empty
if [ -z "$data_folder" ] || [ -z "$dev_pct" ] || [ -z "$max_dev_samples" ] || [ -z "$shuf_seed" ]; then
    echo "Error: All arguments (-i, -p, -m, -s) are required and cannot be empty."
    echo "Usage: $0 -i input_folder -p output_folder -m max_dev_samples -s seed"
    exit 1
fi


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
	# dev_sample_txt=${data_folder}/${basename}.samples.dev.txt
	# train_sample_txt=${data_folder}/${basename}.samples.train.txt

	# Get training/dev beds/txts for each sample bed
	echo "Partitioning samples"
	./code/partition_samples.sh $sample_bed $dev_pct $max_dev_samples $shuf_seed
done
sh ./code/aggregate_and_track_samples.sh $data_folder $shuf_seed 