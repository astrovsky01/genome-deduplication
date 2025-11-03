import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input directory with full kmer counts files", required=True)
args=parser.parse_args()

input_dir=args.input

source_kmer_counts = os.path.join(input_dir, "source_fasta_kmer_counts.txt")
deduped_kmer_counts = os.path.join(input_dir, "deduped_kmer_counts.txt")
ignored_kmer_counts = os.path.join(input_dir, "ignored_kmer_counts.txt")
ignored_and_deduped = os.path.join(input_dir, "ignored_deduped_kmer_counts.txt")
masked_kmer_counts = os.path.join(input_dir, "masked_kmer_counts.txt")
combined_counts = os.path.join(input_dir, "combined_all_kmers.txt")

def read_kmer_counts(file_path):
    try:
        kmer_dict={}
        with open(file_path, 'r') as f:
            for line in f:
                parts=line.strip().split('\t')
                kmer=parts[0]
                count=int(parts[1])
                kmer_dict[kmer]=count
        return kmer_dict
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        sys.exit(1)

source_kmers = read_kmer_counts(source_kmer_counts)
deduped_kmers = read_kmer_counts(deduped_kmer_counts)
ignored_kmers = read_kmer_counts(ignored_kmer_counts)
ignored_and_deduped_kmers = read_kmer_counts(ignored_and_deduped)

def plot_kmer_distributions_grouped(kmer_dicts, labels, colors):
    """Plot multiple k-mer distributions as grouped bar chart"""
    max_count = max(max(d.values()) for d in kmer_dicts)
    
    bin_size = 1
    bins = np.arange(1, max_count + bin_size + 1, bin_size)
    
    all_histograms = []
    for kmer_dict in kmer_dicts:
        counts = list(kmer_dict.values())
        hist, _ = np.histogram(counts, bins=bins)
        all_histograms.append(hist)
    
    num_bins = len(bins) - 1
    num_datasets = len(kmer_dicts)
    bar_width = 0.15
    
    x = np.arange(num_bins)
    for bin_idx in range(num_bins):
        non_zero_bars = []
        for dataset_idx, hist in enumerate(all_histograms):
            if hist[bin_idx] > 0:
                non_zero_bars.append((dataset_idx, hist[bin_idx]))
        
        num_non_zero = len(non_zero_bars)
        if num_non_zero > 0:
            total_width = num_non_zero * bar_width
            start_offset = -total_width / 2
            
            for i, (dataset_idx, height) in enumerate(non_zero_bars):
                offset = start_offset + (i + 0.5) * bar_width
                plt.bar(x[bin_idx] + offset, height, bar_width, 
                       color=colors[dataset_idx], alpha=0.7,
                       edgecolor='white', linewidth=0.5,
                       label=labels[dataset_idx] if bin_idx == 0 else "")
    
    plt.yscale('log')
    plt.xlabel('K-mer Occurence Count')   
    plt.ylabel('Frequency')
    step = max(1, num_bins//20)
    tick_positions = x[::step]
    tick_labels = bins[:-1][::step].astype(int)
    min_len = min(len(tick_positions), len(tick_labels))
    plt.xticks(tick_positions[:min_len], tick_labels[:min_len])
    
    handles = [plt.Rectangle((0,0),1,1, color=colors[i], alpha=0.7) for i in range(num_datasets)]
    plt.legend(handles, labels)
plt.figure(figsize=(12, 6))

source_kmers = read_kmer_counts(source_kmer_counts)
deduped_kmers = read_kmer_counts(deduped_kmer_counts)
ignored_kmers = read_kmer_counts(ignored_kmer_counts)
ignored_and_deduped_kmers = read_kmer_counts(ignored_and_deduped)
masked_kmers = read_kmer_counts(masked_kmer_counts)
combined_kmers = read_kmer_counts(combined_counts)

kmer_dicts = [source_kmers, combined_kmers, ignored_and_deduped_kmers, deduped_kmers, ignored_kmers, masked_kmers]
labels = ["Source", "Combined", "Ignored+Deduped", "Deduped", "Ignored", "Masked"]
colors = ['blue', 'brown', 'orange', 'red', 'green', 'purple']

plot_kmer_distributions_grouped(kmer_dicts, labels, colors)

plt.title('K-mer Count Distributions (Grouped)')
plt.tight_layout()
plt.savefig(os.path.join(input_dir, "kmer_count_distributions.png"))