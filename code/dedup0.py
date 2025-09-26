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

import random

# Convert kmer to int
# Only works up to 32-mers (can be encoded in 64 bits if only A/C/G/T bases)
def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num

# Not necessarily used, but good for testing
def decode_kmer(kmer_num, k=32):
    char_map = {0:'A', 1:'C', 2:'G', 3:'T'}
    kmer = []
    for i in range(k):
        nucleotide_code = kmer_num & 3 # 0x11
        kmer.append(char_map[nucleotide_code])
        kmer_num >>= 2
    return ''.join(reversed(kmer))

# sample_regions are inclusive at the start, exclusive at the end
# masked_starts are the start coordinates of the masked kmers
# skipped_regions are inclusive at the start, exclusive at the end
def deduplicate(seq, k, sample_len, allowed_duplicate_rate=0.0, seen_kmers=None, min_sample_len=None, seed=123):

    rng = random.Random(seed)

    # Check validity of inputs
    if k > sample_len:
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
        gap_index = seq[sample_start:sample_end].find('N')        
        if gap_index > -1:
            sample_end = sample_start + gap_index + 1
            skipped_regions.append((sample_start, sample_end))
            sample_start = sample_end
            continue

        # Loop through all kmers in this possible sample
        for kmer_start_idx in range(sample_start, sample_end-k+1):

            # Get kmer and its encoding for storage
            kmer = seq[kmer_start_idx:kmer_start_idx+k]
            kmer_num = encode_kmer(kmer)

            # If we haven't seen this kmer before, record it in the sample kmers
            if kmer_num not in seen_kmers and kmer_num not in sample_seen_kmers:
                sample_seen_kmers.add(kmer_num)

            # If we've seen this kmer before, stop analyzing this sample
            else:

                # Stochasticity is built in here here; if random() < allowed duplicate rate, let through this duplicate
                if rng.random() < allowed_duplicate_rate:
                    continue
                
                # Handle finding a repeat
                masked_starts.append(kmer_start_idx)
                sample_end = kmer_start_idx
                exited_early = True
                break
       
        # Check if we found a valid sample (didn't exit too early)
        # If exited early, sample_end is sample_end + k - 1
        sample_end_idx = sample_end - 1 + k if exited_early else sample_end
        if (sample_end_idx - sample_start) >= min_sample_len:

            # Add all sample seen kmers to seen kmers
            seen_kmers.update(sample_seen_kmers) 

            # Record sample indices
            sample_regions.append((sample_start, sample_end_idx))

            # TODO: add ignored kmers at the edges of samples?

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

# Test encoder/decoder
kmers = ["GATTACA", "GATTACAT", "TCCATGGAC"]
for kmer in kmers:
    print(f"kmer = {kmer}; encoded kmer = {encode_kmer(kmer)}; decoded kmer = {decode_kmer(encode_kmer(kmer), len(kmer))}")
exit()

#sequence = "AAACCCAACACCGGGGGGTGTGTGAAA"
#sequence = "AAACCCAACACC"
#sequence = "AAACCCAAACCC"
sequence = "AAANAAACCCAACACCNGGGGGNT"
k = 3
sample_len = 4

sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate(sequence, k, sample_len)

print(f"Sample regions: {sample_regions}")
print(f"Masked starts: {masked_starts}")
print(f"Skipped regions: {skipped_regions}")
print(f"Seen kmers: {seen_kmers}")

