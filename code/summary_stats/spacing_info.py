import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math

import argparse

parser = argparse.ArgumentParser(description="Calculate spacing information between deduplicated regions and plot histograms.")
parser.add_argument("--input_file", default="output_", help="Input csv file")
parser.add_argument("--type", default="dedup", help="Type of spacing information", choices=["dedup", "ignored", "masked"])
parser.add_argument("--output_dir", default=".", help="Directory to save output plots")
parser.add_argument("--bincount", type=int, default=100, help="Number of bins for the histogram")
args = parser.parse_args()

def plot_lengths(data, outname, title, bins=20):
    plt.figure(figsize=(10, 6))
    counts, bin_edges, patches = plt.hist(data, bins=bins, edgecolor='black')
    count_min = np.min(counts)  # Minimum non-zero count
    count_max = np.max(counts)
    count_range = math.log(count_max - count_min)
    
    # Use log scale if the range spans more than 2 orders of magnitude
    if count_range > 2:
        plt.yscale("log")
        # Set y-limits to show all data including 1-count bars - use 0.1 to ensure 1 is fully visible
        plt.ylim(0.1, count_max * 1.1)
        ylabel = "Count (log)"
    else:
        plt.yscale("linear")
        # Set y-limits starting from 0 for linear scale
        plt.ylim(0, count_max * 1.1)
        ylabel = "Count"
    
    # FORCE plain formatting
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    ax.xaxis.set_minor_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    
    plt.xlabel("Length (bp)")
    plt.ylabel(ylabel)
    plt.title(title)
    
    # Add range information as text on the plot
    x_min = np.min(data)
    x_max = np.max(data)
    range_text = f"Range: {x_min:,} - {x_max:,} bp"
    plt.text(0.02, 0.98, range_text, transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.savefig(outname + ".png")
    plt.close()

if __name__ == "__main__":
    input_file = args.input_file
    output_dir = args.output_dir
    bincount = args.bincount
    basename = os.path.basename(input_file).replace(".csv", "")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    spacing_data = np.loadtxt(input_file, delimiter=',', dtype=int)
    type=args.type
    if args.type == "dedup":
        plot_lengths(spacing_data, os.path.join(output_dir, f"{basename}_distance_between_deduplicated_regions"), 
                     "Distances Between Deduplicated Regions", bins=bincount)
    # Plot histogram of spacing data
    elif args.type == "ignored":
        plot_lengths(spacing_data, os.path.join(output_dir, f"{basename}_size_of_ignored_regions"), 
                     "Lengths of Ignored Regions", bins=bincount)
    elif args.type == "masked":
        plot_lengths(spacing_data, os.path.join(output_dir, f"{basename}_size_of_masked_regions"), 
                     "Lengths of Masked Regions", bins=bincount)