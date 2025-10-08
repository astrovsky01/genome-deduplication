#!/usr/bin/env bash

### Usage: ./code/check_beds.sh [output dir from dedup0.py] [sample prefix (everything before the .samples.bed)]
### e.g. ./code/check_beds.sh tests/out/new_no_overlap/ GCA_001775245.1_ASM177524v1_genomic

# Input args
data_folder=$1
sample=$2

### Check 4 - beds are fully disjoint
bed_types=("samples" "masks" "ignored" "ambiguous")
n_bed_types=${#bed_types[@]}
success=1
for ((i=1; i<$n_bed_types; i++)); do
	for ((j=0; j<$i; j++)); do
		bed_A=$data_folder/${sample}.${bed_types[$i]}.bed
		bed_B=$data_folder/${sample}.${bed_types[$j]}.bed
		n_overlaps=$(bedtools intersect -a $bed_A -b $bed_B | wc -l)
		if [[ $n_overlaps -gt 0 ]]; then
			success=0
        		echo "Check 4 failed: ${bed_types[$i]} and ${bed_types[$j]} bed files have $n_overlaps overlapping regions"
		fi
	done
done
if [[ $success -eq 1 ]]; then
	echo "Check 4 complete! All bed files are disjoint"
fi

### Check 5 - beds cover the entirety of the original fasta

# Get corresponding fasta
sample_fa=$(grep "^$sample" ${data_folder}/basename_fasta_match.txt | cut -f2)

# Get chromosome sizes 
chrom_size_bed=chrom_sizes.bed
if [[ "$sample_fa" == *.gz ]]; then
	zcat $sample_fa | awk 'BEGIN {OFS="\t"} /^>/ {if (len) {print name, '0', len}; name=substr($1,2); len=0; next} {len+=length($0)} END {print name, '0', len}' >> $chrom_size_bed
else
	awk 'BEGIN {OFS="\t"} /^>/ {if (len) {print name, '0', len}; name=substr($1,2); len=0; next} {len+=length($0)} END {print name, '0', len}' $sample_fa >> $chrom_size_bed
fi

# Merge sample, masked, and ignored beds
merged_regions_bed=merged.bed
cat $data_folder/${sample}.*.bed | sortBed | mergeBed > $merged_regions_bed

# Check that the merged regions cover the entirety of all chromosomes
n_diff_lines=$(diff $chrom_size_bed $merged_regions_bed | wc -l)
if [[ $n_diff_lines -eq 0 ]]; then
	echo "Check 5 complete! The outputted regions encompass the entirety of all chromosomes"
else
	echo "Check 5 failed: up to $n_diff_lines chromosomes are not fully encompassed by the outputted regions"
fi

# Clean up by removing intermediate files
rm -f $chrom_size_bed
rm -f $merged_regions_bed
