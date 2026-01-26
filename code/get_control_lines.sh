#!/bin/bash
#Get the actual lines from deduped control files based on control line counts

# Parse command line arguments
while getopts "c:d:s:o:" opt; do
    case $opt in
        c)
            control_dir="$OPTARG"
            ;;
        d)
            deduped_lines_dir="$OPTARG"
            ;;

        s)
            sample="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
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
if [ -z "$control_dir" ] || [ -z "$deduped_lines_dir" ] || [ -z "$sample" ] || [ -z "$output_dir" ]; then
    echo "Error: -c (control directory), -d (deduped lines directory), -s (sample), and -o (output directory) are required."
    echo "Usage: $0 -c <control_directory> -d <deduped_lines_directory> -s <sample> -o <output_directory>"
    exit 1
fi

# Check if directories exist
if [ ! -d "$control_dir" ]; then
    echo "Error: Control directory does not exist: $control_dir"
    exit 1
fi

if [ ! -d "$deduped_lines_dir" ]; then
    echo "Error: Deduped lines directory does not exist: $deduped_lines_dir"
    exit 1
fi
train_lines=(${deduped_lines_dir}/${output_dir}/${sample}_train_sampled.lines)
dev_lines=(${deduped_lines_dir}/${output_dir}/${sample}_dev_sampled.lines)

# dev lines
# echo "Generating control dev bed for sample: $sample"
# awk 'NR==FNR{lines[NR]=$1; count=NR; next} {data[FNR]=$0} END{for(i=1; i<=count; i++) print data[lines[i]]}' "$dev_lines" "${control_dir}/all_dev.bed" > "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.bed"
# dev_output_line_count=$(wc -l < "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.bed")
# dev_input_line_count=$(wc -l < "${deduped_lines_dir}/all_dev.bed")
# if [ "$dev_input_line_count" -eq "$dev_output_line_count" ]; then
#     echo "✓ Dev files match: $dev_input_line_count lines"
# else
#     echo "✗ Dev line count mismatch: input=$dev_input_line_count, output=$dev_output_line_count"
#     exit 1
# fi

# echo "Generating control dev txt for sample: $sample"
# awk 'NR==FNR{lines[NR]=$1; count=NR; next} {data[FNR]=$0} END{for(i=1; i<=count; i++) print data[lines[i]]}' "$dev_lines" "${control_dir}/all_dev.txt" > "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.txt"
# dev_output_line_count=$(wc -l < "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.txt")
# dev_input_line_count=$(wc -l < "${deduped_lines_dir}/all_dev.txt")
# if [ "$dev_input_line_count" -eq "$dev_output_line_count" ]; then
#     echo "✓ Dev files match: $dev_input_line_count lines"
# else
#     echo "✗ Dev line count mismatch: input=$dev_input_line_count, output=$dev_output_line_count"
#     exit 1
# fi

# # train lines
# echo "Generating control train bed for sample: $sample"
# awk 'NR==FNR{lines[NR]=$1; count=NR; next} {data[FNR]=$0} END{for(i=1; i<=count; i++) print data[lines[i]]}' "$train_lines" "${control_dir}/all_train.bed" > "${deduped_lines_dir}/${output_dir}/${sample}_control_train.bed"
# train_output_line_count=$(wc -l < "${deduped_lines_dir}/${output_dir}/${sample}_control_train.bed")
# train_input_line_count=$(wc -l < "${deduped_lines_dir}/all_train.bed")
# if [ "$train_input_line_count" -eq "$train_output_line_count" ]; then
#     echo "✓ Train files match: $train_input_line_count lines"
# else
#     echo "✗ Train line count mismatch: input=$train_input_line_count, output=$train_output_line_count"
#     exit 1
# fi
# echo "Generating control train txt for sample: $sample"
# awk 'NR==FNR{lines[NR]=$1; count=NR; next} {data[FNR]=$0} END{for(i=1; i<=count; i++) print data[lines[i]]}' "$train_lines" "${control_dir}/all_train.txt" > "${deduped_lines_dir}/${output_dir}/${sample}_control_train.txt"
# train_output_line_count=$(wc -l < "${deduped_lines_dir}/${output_dir}/${sample}_control_train.txt")
# train_input_line_count=$(wc -l < "${deduped_lines_dir}/all_train.txt")
# if [ "$train_input_line_count" -eq "$train_output_line_count" ]; then
#     echo "✓ Train files match: $train_input_line_count lines"
# else
#     echo "✗ Train line count mismatch: input=$train_input_line_count, output=$train_output_line_count"
#     exit 1
# fi  
# echo "Control line extraction completed successfully for sample: $sample"

# Memory-efficient extraction for terabyte files
# This approach uses temporary files to avoid loading entire source file into memory

echo "Starting memory-efficient extraction for terabyte files..."

temp_dir="/tmp/control_extract_$$"
mkdir -p "$temp_dir"

# Dev files - memory efficient approach  
echo "Processing dev files with memory-efficient method..."

# Create temp file with original positions and line numbers
nl -nln "$dev_lines" > "$temp_dir/dev_with_positions.txt"

# Sort by line number (column 2)
sort -k2,2n "$temp_dir/dev_with_positions.txt" > "$temp_dir/dev_sorted.txt"

# Extract in sorted order (single pass through huge file)
awk '{print $2}' "$temp_dir/dev_sorted.txt" | awk 'NR==FNR{a[$1]=1;next} FNR in a' - "${control_dir}/all_dev.bed" > "$temp_dir/dev_extracted.bed"
awk '{print $2}' "$temp_dir/dev_sorted.txt" | awk 'NR==FNR{a[$1]=1;next} FNR in a' - "${control_dir}/all_dev.txt" > "$temp_dir/dev_extracted.txt"

# Re-sort back to original order using positions
paste "$temp_dir/dev_sorted.txt" "$temp_dir/dev_extracted.bed" | sort -k1,1n | cut -f3- > "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.bed"
paste "$temp_dir/dev_sorted.txt" "$temp_dir/dev_extracted.txt" | sort -k1,1n | cut -f3- > "${deduped_lines_dir}/${output_dir}/${sample}_control_dev.txt"

# Train files - memory efficient approach
echo "Processing train files with memory-efficient method..."

# Create temp file with original positions and line numbers
nl -nln "$train_lines" > "$temp_dir/train_with_positions.txt"

# Sort by line number (column 2) 
sort -k2,2n "$temp_dir/train_with_positions.txt" > "$temp_dir/train_sorted.txt"

# Extract in sorted order (single pass through huge file)
awk '{print $2}' "$temp_dir/train_sorted.txt" | awk 'NR==FNR{a[$1]=1;next} FNR in a' - "${control_dir}/all_train.bed" > "$temp_dir/train_extracted.bed"
awk '{print $2}' "$temp_dir/train_sorted.txt" | awk 'NR==FNR{a[$1]=1;next} FNR in a' - "${control_dir}/all_train.txt" > "$temp_dir/train_extracted.txt"

# Re-sort back to original order using positions
paste "$temp_dir/train_sorted.txt" "$temp_dir/train_extracted.bed" | sort -k1,1n | cut -f3- > "${deduped_lines_dir}/${output_dir}/${sample}_control_train.bed"
paste "$temp_dir/train_sorted.txt" "$temp_dir/train_extracted.txt" | sort -k1,1n | cut -f3- > "${deduped_lines_dir}/${output_dir}/${sample}_control_train.txt"

# Clean up temp files
rm -rf "$temp_dir"

echo "Memory-efficient extraction completed successfully for sample: $sample"