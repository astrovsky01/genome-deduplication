#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input arg
data_folder=$1

# Global args
dev_pct=0.1
max_dev_samples=100000
shuf_seed=123
no_txt=false
test_pct=0


while getopts ":i:p:m:s:r:t" opt; do
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
        r)
            test_pct="$OPTARG"
            ;;
        t)
            no_txt=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i input_folder -p dev_pct -m max_dev_samples -s seed [-t]"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo "Usage: $0 -i input_folder -p dev_pct -m max_dev_samples -s seed [-t]"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided and not empty
if [ -z "$data_folder" ] || [ -z "$dev_pct" ] || [ -z "$max_dev_samples" ] || [ -z "$shuf_seed" ]; then
    echo "Error: Required arguments (-i, -p, -m, -s) cannot be empty."
    echo "Usage: $0 -i input_folder -p output_folder -m max_dev_samples -s seed [-t]"
    exit 1
fi

if ! [[ "$test_pct" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]] ||
   ! awk -v pct="$test_pct" 'BEGIN { exit !(pct >= 0 && pct <= 1) }'; then
    echo "Error: Test percentage must be a number between 0 and 1." >&2
    exit 1
fi

if ! [[ "$max_dev_samples" =~ ^[0-9]+$ ]]; then
    echo "Error: Maximum dev samples must be a non-negative integer." >&2
    exit 1
fi

if ! [[ "$shuf_seed" =~ ^-?[0-9]+$ ]]; then
    echo "Error: Shuffle seed must be an integer." >&2
    exit 1
fi

if ! [[ "$dev_pct" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]] ||
   ! awk -v pct="$dev_pct" 'BEGIN { exit !(pct >= 0 && pct <= 1) }'; then
    echo "Error: Dev percentage must be a number between 0 and 1." >&2
    exit 1
fi


#Check if paritioned files exist already and check for confirmation to overwrite
if ls ${data_folder}/*.samples.dev.bed 1> /dev/null 2>&1 || ls ${data_folder}/*.samples.train.bed 1> /dev/null 2>&1; then
    read -p "Partitioned sample files already exist. Do you want to overwrite them? (y/n) " choice
    case "$choice" in
        y|Y ) echo "Overwriting existing partitioned sample files..." && rm -f "${data_folder}"/*.samples.dev.bed "${data_folder}"/*.samples.train.bed "${data_folder}"/*.samples.test.bed;;
        n|N ) echo "Exiting script."; exit 1;;
        * ) echo "Invalid choice. Exiting script."; exit 1;;
    esac
fi

#Check if final output files exist already and check for confirmation to overwrite
if [ -f "${data_folder}/all_train.bed" ] || [ -f "${data_folder}/all_dev.bed" ] || [ -f "${data_folder}/all_test.bed" ]; then
    read -p "Combined and shuffled files already exist. Do you want to overwrite them? (y/n) " choice
    case "$choice" in
        y|Y ) echo "Overwriting existing files..." && rm -f "${data_folder}/all_train.bed" "${data_folder}/all_dev.bed" "${data_folder}/all_test.bed";;
        n|N ) echo "Exiting script."; exit 1;;
        * ) echo "Invalid choice. Exiting script."; exit 1;;
    esac
fi

if [ -f "${data_folder}/all_train.txt" ] || [ -f "${data_folder}/all_dev.txt" ] || [ -f "${data_folder}/all_test.txt" ]; then
    read -p "Combined and shuffled text files already exist. Do you want to overwrite them? (y/n) " choice
    case "$choice" in
        y|Y ) echo "Overwriting existing text files..." && rm -f "${data_folder}/all_train.txt" "${data_folder}/all_dev.txt" "${data_folder}/all_test.txt";;
        n|N ) echo "Exiting script."; exit 1;;
        * ) echo "Invalid choice. Exiting script."; exit 1;;
    esac
fi



## Create train/dev txt files for each sample in the input directory

# Get all samples.bed files in the input directory
sample_beds=$(ls $data_folder/*.samples.bed)

for sample_bed in $sample_beds; do
    if [[ "$(basename "$sample_bed")" == "shuffled.samples.bed" ]]; then
        continue
    fi

    # Get all needed file names
    basename=$(basename "$sample_bed" | rev | cut -d'.' -f3- | rev)
    # sample_fasta=$(grep "^${basename}" "${data_folder}/basename_fasta_match.txt" | cut -f2)
    # dev_sample_bed="${data_folder}/${basename}.samples.dev.bed"
    # train_sample_bed="${data_folder}/${basename}.samples.train.bed"
    # test_sample_bed="${data_folder}/${basename}.samples.test.bed"

    echo "Partitioning samples"
    "$SCRIPT_DIR/partition_samples.sh" "$sample_bed" "$dev_pct" "$max_dev_samples" "$shuf_seed" "$test_pct"
done
if [ "$no_txt" = false ]; then
    echo "Finished partitioning samples into train/dev sets. Aggregating and retrieving sequences now"
    "$SCRIPT_DIR/aggregate_and_track_samples.sh" "$data_folder" "$shuf_seed" "$test_pct"
else
    echo "Finished paritioning samples into train/dev sets. Aggregating and shuffling bed files."
    "$SCRIPT_DIR/aggregate_and_track_samples.sh" "$data_folder" "$shuf_seed" "$test_pct" true
fi
