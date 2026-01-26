#!/bin/bash

# Usage: ./generate_controls.sh -c <control_directory> -d <deduped_directory>

# Parse command line arguments
while getopts "c:d:" opt; do
    case $opt in
        c)
            control_dir="$OPTARG"
            ;;
        d)
            deduped_dir="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

# Check if both arguments are provided
if [ -z "$control_dir" ] || [ -z "$deduped_dir" ]; then
    echo "Error: Both -c (control directory) and -d (deduped directory) are required."
    echo "Usage: $0 -c <control_directory> -d <deduped_directory>"
    exit 1
fi

# Check if directories exist
if [ ! -d "$control_dir" ]; then
    echo "Error: Control directory does not exist: $control_dir"
    exit 1
fi

if [ ! -d "$deduped_dir" ]; then
    echo "Error: Deduped directory does not exist: $deduped_dir"
    exit 1
fi


#Find line counts for control files
control_train=("$control_dir"/all_train.bed)
control_dev=("$control_dir"/all_dev.bed)
control_train_lines=$(wc -l < "${control_train[0]}")
control_dev_lines=$(wc -l < "${control_dev[0]}")
echo "Control train lines: $control_train_lines"
echo "Control dev lines: $control_dev_lines"
sample_name=$(basename "$deduped_dir")
#Generate deduped control line number files
python get_samples_based_on_dedup.py -f "$deduped_dir" -t "$control_train_lines" -d "$control_dev_lines" -s 42 -o "control_samples"
sh get_control_lines.sh -c "$control_dir" -d $deduped_dir -o "control_samples" -s "$sample_name"