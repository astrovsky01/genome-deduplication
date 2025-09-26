### New attempt to deduplicate a fasta
### Deduplicates and samples at the same time
### Only consider a kmer if it is within a sample
### Basic algorithm:
###     1. Sample sequence to first `sample_len` bases
###     2. Iterate over all kmers in sample
###     3. If kmer not in seen_kmers, add it to seen_kmers
###     4. If kmer in seen_kmers, jump to next possible sample and repeat
### Notes:
###     1. Most issues arise when we find a seen kmer in the sample. If we do,
###        we need to do a few things. For one, we jump to the next possible
###        sample, which starts at the current index + 1. Doing so should
###        retroactively invalidate all kmers we saw in that failed sample.
###        To implement this, we shouldn't add the sample's kmers to the seen
###        kmer set until we finalize the sample
### This program basically divides a fasta into three categories:
###     1. Samples - valid regions with no repeats, on which we train
###     2. Masked regions - regions containing duplicates
###     3. Ignored regions - regions that we threw out because we found a
###        duplicate later on in the sample region
###        e.g. seq = AAATTAAATTAA, k=3, sample_len=4
###        Consider AAAT sample
###        Look at AAA, add to seen_kmers
###        Look at AAT, add to seen_kmers
###        No repeats in sample, so add AATA to samples
###        Consider TAAA sample
###        Look at TAA, add to seen_kmers
###        Look at AAA; it's seen before
###        Add AAA to masked regions
###        Add TAA to ignored regions
###        Jump ahead to next possible sample, TTAA
### The entirety of the input should fall into 1 of the 3 categories

# sample_regions are inclusive at the start, exclusive at the end
# masked_starts are the start coordinates of the masked kmers
# skipped_regions are inclusive at the start, exclusive at the end

import argparse
import Bio
from Bio import SeqIO
import pickle
import os
import random


## Variable setting ================
parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="Input FASTA file")
parser.add_argument("-k". "--kmer", type=int, default=32, help="Kmer size (default: 32)")
parser.add_argument("-s", "--sample_len", type=int, default=1000, help="Sample length (default: 1000)")
parser.add_argument("-m", "--min_sample_len", type=int, default=50, help="Minimum sample length (default: 50)")
parser.add_argument("-p", "--seen_kmers", default=None, help="Pickle file containing seen kmers (default: None)")
parser.add_argument("-r", "--retain", type=float, default=0.0, help="Likelihood a duplicate kmer will be allowed through as a value from [0,1]")

args = parser.parse_args()

fasta = args.fasta
file_basename = os.path.splitext(os.path.basename(fasta))[0]
k = args.k
sample_len = args.sample_len
min_sample_len = args.min_sample_len
if args.seen_kmers is not None:
    seen_kmers = pickle.load(open(args.seen_kmers, "rb"))
else:
    seen_kmers = None
if args.retain is not None:
    if args.stochastic_retention > 1.0 or args.retain < 0.0:
        raise("Error: Retention rate of duplicate kmers cannot be greater than 1.0 or less than 0.0")
    else:
        retain = args.retain
local_dict = {}

# =====================================


## Component functions ================

def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num


def n_check(seq, sample_start, sample_end, skipped_regions):
    "Skip over regions with N"
    gap_index = seq[sample_start:sample_end].find('N')        
    if gap_index > -1:
        sample_end = sample_start + gap_index + 1
        skipped_regions.append((sample_start, sample_end))
        sample_start = sample_end
    return(skipped_regions, sample_start, sample_end)

def sample_scan(seq, start, end, k, seen_kmers, sample_seen_kmers, masked_starts, dedup_retain=globals()["retain"]):
    """Scan a sample for kmers. If a duplicate is found either within the same set or globally, 
    return the index of the start of the duplicate kmer. If no duplicate is found, """
    for kmer_start_idx in range(start, end-k+1):

        # Get encoded kmer
        kmer = encode_kmer(seq[kmer_start_idx:kmer_start_idx+k])

        # If we haven't seen this kmer in this sample before, record it in the sample kmers
        if kmer not in seen_kmers and kmer not in sample_seen_kmers:
            sample_seen_kmers.add(kmer)
        # If we've seen this kmer before, stop analyzing this sample unless it passes random check
        else:
            if random.random() > dedup_retain:
                continue
            else:
                # Handle finding a repeat
                return sample_seen_kmers, masked_starts, True, kmer_start_idx
    return sample_seen_kmers, masked_starts, False, end

# =====================================


## Output Functions ================

def writeout_bed(name, regions, bedfile):
    """
    Create bed files for regions or just a list of start sites for beginning positions
    """
    with open(bedfile, 'w') as f:
        if len(regions[0]) == 2:
            for start, end in regions:
                f.write(f"{name}\t{start}\t{end}\n")
        else:
            for start in regions:
                f.write(f"{name}\t{start}\n")

def writeout_kmers(seen_kmer_set, outfile):
    """
    Dump seen kmers into a pickle file that can be used as input for next dedup file
    """
    with open(outfile, 'w') as f:
        pickle.dump(seen_kmer_set, f)

def output_dump(local_dict, file_basename):
    """
    Output all bed files
    """
    with open(file_basename + "_sample_regions.bed", 'w') as sample_regions:
        with open(file_basename + "_masked_starts.bed", 'w') as masked_starts:
            with open(file_basename + "_skipped_regions.bed", 'w') as skipped_regions:
                for seqname in local_dict.keys():
                    writeout_bed(seqname, local_dict[seqname]["sample_regions"], sample_regions)
                    writeout_bed(seqname, local_dict[seqname]["masked_starts"], masked_starts)
                    writeout_bed(seqname, local_dict[seqname]["skipped_regions"], skipped_regions)
# =====================================


## Main Deduplication ================
def deduplicate(seq, k, sample_len, seen_kmers=None, min_sample_len=None):

    # Check validity of inputs
    if k < sample_len:
        raise("Error: sample length cannot be greater than k.")
    if min_sample_len is not None and len(seq) < min_sample_len:
        raise("Error: the sequence length is less than min_sample_len.")
    if min_sample_len is None and len(seq) < sample_len:
        raise("Error: the sequence length is less than sample_len.")

    # Instantiate seen_kmers if None
    if seen_kmers is None:
        seen_kmers = set()

    # Set min_sample_len to sample_len if None
    if min_sample_len is None:
        min_sample_len = sample_len

    # Set up storage for samples, masked regions, and skipped regions
    sample_regions = []
    masked_starts = []
    skipped_regions = []

    # Set up indices of all possible samples 
    max_start_idx = len(seq) - sample_len
    sample_start = 0

    # Investigate every possible sample
    while sample_start <= max_start_idx:

        # Get boundary for this possible sample
        sample_end = sample_start + sample_len

        # Data structure for new kmers in this sample - see note above
        sample_seen_kmers = set()

        # Record whether we exited early (found a repeat kmer)
        exited_early = False

        # Early exit if N is in this possible sample
        # Although code-wise, don't say exited_early
        # gap_index = seq[sample_start:sample_end].find('N')        
        # if gap_index > -1:
        #     sample_end = sample_start + gap_index + 1
        #     skipped_regions.append((sample_start, sample_end))
        #     sample_start = sample_end
        #     continue
        skipped_regions, sample_start, sample_end = n_check(seq, sample_start, sample_end, skipped_regions)

        # Loop through all kmers in this possible sample
        # for kmer_start_idx in range(sample_start, sample_end-k+1):

        #     # Get kmer
        #     kmer = seq[kmer_start_idx:kmer_start_idx+k]

        #     # If we haven't seen this kmer before, record it in the sample kmers
        #     if kmer not in seen_kmers and kmer not in sample_seen_kmers:
        #         sample_seen_kmers.add(kmer) # TODO: should encode kmers

        #     # If we've seen this kmer before, stop analyzing this sample
        #     else:

        #         # TODO: Could add stochasticity here; if random() < dedup_pct; continue
                
        #         # Handle finding a repeat
        #         masked_starts.append(kmer_start_idx)
        #         sample_end = kmer_start_idx
        #         exited_early = True
        #         break
        sample_seen_kmers, masked_starts, exited_early, sample_end = sample_scan(seq, sample_start, sample_end, k, seen_kmers, sample_seen_kmers, masked_starts)
       
        # Check if we found a valid sample (didn't exit too early)
        # If exited early, sample_end is sample_end + k - 1
        sample_end_idx = sample_end - 1 + k if exited_early else sample_end
        if (sample_end_idx - sample_start) >= min_sample_len:

            # Add all sample seen kmers to seen kmers
            seen_kmers.update(sample_seen_kmers) 

            # Record sample indices
            sample_regions.append((sample_start, sample_end_idx))

            # TODO: add ignored kmers at the edges of samples

        # If not a valid sample, record skipped region
        else:

            # Write skipped region
            # If sample is too short, we want to use sample_end, not sample_end_idx
            skipped_regions.append((sample_start, sample_end))
    
        # Update sample_start to look at next possible sample            
        if exited_early: # If so, sample_end is 1 short (start of repeat kmer)
            sample_end += 1
        sample_start = sample_end

    # Handle possible skipped sequence at the end
    # If sample_start is more than max_start_idx, we skipped ahead at some point
    # Therefore we did not consider the final possible sample starting at max_start_idx
    if not exited_early and sample_start > max_start_idx and sample_start < len(seq):
        skipped_regions.append((sample_start, len(seq)))

    return sample_regions, masked_starts, skipped_regions, seen_kmers
# =====================================




# #sequence = "AAACCCAACACCGGGGGGTGTGTGAAA"
# #sequence = "AAACCCAACACC"
# #sequence = "AAACCCAAACCC"
# sequence = "AAANAAACCCAACACCNGGGGGNT"
# k = 3
# sample_len = 4






with open(fasta, 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        sequence = str(record.seq)
        seqname = record.id
        sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate(sequence, k, sample_len, seen_kmers=seen_kmers)

        print(f"Sample regions: {sample_regions}")
        print(f"Masked starts: {masked_starts}")
        print(f"Skipped regions: {skipped_regions}")

        # Keep dict associating the regions with the sequence name
        local_dict[seqname] = {
            "sample_regions": sample_regions,
            "masked_starts": masked_starts,
            "skipped_regions": skipped_regions,
        }
        print(f"Seen kmers: {seen_kmers}")
    output_dump(local_dict, file_basename)
writeout_kmers(seen_kmers, file_basename + f".pickle")

