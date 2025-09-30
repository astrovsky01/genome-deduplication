## Component functions ================

import argparse
import Bio
from Bio import SeqIO
import pickle
import os
import random

def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num

def decode_kmer(kmer_num, k):
    char_map = {0:'A', 1:'C', 2:'G', 3:'T'}
    kmer = []
    for i in range(k):
        nucleotide_code = kmer_num & 3 # 0x11
        kmer.append(char_map[nucleotide_code])
        kmer_num >>= 2
    return ''.join(reversed(kmer))

def n_check(seq, sample_start, sample_end, skipped_regions):
    "Skip over regions with N"
    
    gap_index = seq[sample_start:sample_end][::-1].find('N')        
    if gap_index > -1:
        gap_index = sample_end - (sample_start + gap_index) - 1
        sample_end = sample_start + gap_index + 1
        print(gap_index)
        skipped_regions.append((sample_start, sample_end))
        sample_start = sample_end
    return(skipped_regions, sample_start, sample_end)

def sample_scan(seq, start, end, k, seen_kmers, sample_seen_kmers, masked_starts, dedup_retain=0.0, rng=None):
    """Scan a sample for kmers. If a duplicate is found either within the same set or globally, 
    return the index of the start of the duplicate kmer. If no duplicate is found, """
    if rng is None:
        rng = random.Random(123)
        
    for kmer_start_idx in range(start, end-k+1):

        # Get encoded kmer
        kmer = encode_kmer(seq[kmer_start_idx:kmer_start_idx+k])

        # If we haven't seen this kmer in this sample before, record it in the sample kmers
        if kmer not in seen_kmers and kmer not in sample_seen_kmers:
            sample_seen_kmers.add(kmer)
        # If we've seen this kmer before, stop analyzing this sample unless it passes random check
        else:
            if rng.random() > dedup_retain:
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
        # Region file vs start site file
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
    with open(outfile, 'wb') as f:
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