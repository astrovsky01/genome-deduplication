#!/usr/bin/env bash

# Usage: ./code/get_sample_seqs.sh [fasta] [beds...]
# e.g.:  ./code/get_sample_seqs.sh human.fa human.dev.bed human.train.bed
# Produces a corresponding txt file for each input bed
# e.g. for the above example, produces human.dev.txt and human.train.txt

module load bedtools
module load samtools

#in_bed=tests/out/dummy/dummy_1.samples.bed
#in_fa=tests/test-data/dummy_1.fa
#out_file=foo.txt

#in_bed=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.bed
#in_fa=/scratch4/mschatz1/mrevsin1/ncbi_genomes/genomes/GCA_001775245.1_ASM177524v1_genomic.fna.gz
#out_file=tests/out_small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.txt

in_fa=$1
in_beds=${@:2}

in_fa_ext="${in_fa##*.}"
if [[ "$in_fa_ext" == "gz" ]]; then
	zcat $in_fa > foo.fa
	samtools faidx foo.fa
	for in_bed in ${in_beds[@]}; do
		out_file=${in_bed::-4}.txt
		bedtools getfasta -fi foo.fa -bed $in_bed -tab | cut -f2 > $out_file
	done
	rm -f foo.fa foo.fa.fai
else
	for in_bed in ${in_beds[@]}; do
		out_file=${in_bed::-4}.txt
		bedtools getfasta -fi foo.fa -bed $in_bed -tab | cut -f2 > $out_file
	done
fi
