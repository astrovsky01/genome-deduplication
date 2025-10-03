# Implement all checks in tests/checks.txt
# While this doesn't guarantee "correct" deduplication, it does cross off potential sources of error
# Usage: python check_validity.py [directory with train/dev txt/bed] [final seen kmer pickle file] [...k]

###============================================================================
### Imports

import os
import pickle
import sys

###============================================================================
### Functions

def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num

###============================================================================
### Check all input args

if len(sys.argv) < 3:
    raise("Error: not enough arguments. Usage: python check_validity.py [data dir] [seen kmer file]")

data_dir = sys.argv[1]
seen_kmer_file = sys.argv[2]
k = 32 if len(sys.argv) == 3 else int(sys.argv[3])

if not os.path.isdir(data_dir):
    print(f"Error: could not find data directory at {data_dir}")
    exit()

if not os.path.isfile(seen_kmer_file):
    print(f"Error: could not find seen kmer file at {seen_kmer_file}")
    exit()

train_txt = os.path.join(data_dir, "train.txt")
train_bed = os.path.join(data_dir, "train.bed")
dev_txt = os.path.join(data_dir, "dev.txt")
dev_bed = os.path.join(data_dir, "dev.bed")

if not (os.path.isfile(train_txt) and os.path.isfile(train_bed) and os.path.isfile(dev_txt) and os.path.isfile(dev_bed)):
    raise("Error: one or more of train.txt, train.bed, dev.txt, or dev.bed is not found in the supplied data directory")

###============================================================================
### Read in all data

with open(train_txt, 'r') as f:
    train_samples = [line.rstrip('\n') for line in f.readlines()]
with open(train_bed, 'r') as f:
    train_coords = [line.rstrip('\n').split('\t') for line in f.readlines()]
with open(dev_txt, 'r') as f:
    dev_samples = [line.rstrip('\n') for line in f.readlines()]
with open(dev_bed, 'r') as f:
    dev_coords = [line.rstrip('\n').split('\t') for line in f.readlines()]
with open(seen_kmer_file, 'rb') as f:
    seen_kmers = pickle.load(f)

samples = train_samples + dev_samples

###============================================================================
### Check 1: all kmers of all samples are unique

sample_kmers = [seq[i:i+k] for seq in samples for i in range(len(seq)-k+1)]
assert(len(sample_kmers) == len(set(sample_kmers)))
print("Check 1 complete! All kmers of all samples are unique")

###============================================================================
### Check 2: all kmers of all samples is the same as seen kmers from file

sample_kmer_nums = {encode_kmer(kmer) for kmer in sample_kmers}

assert(len(sample_kmer_nums.difference(seen_kmers)) == 0 and len(seen_kmers.difference(sample_kmer_nums)) == 0)
print("Check 2 complete! All kmers of all samples is the same set as seen kmers")

###============================================================================
### Check 3: all kmers of all samples is the same as seen kmers from file





