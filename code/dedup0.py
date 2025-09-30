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

###=============================================================================
### Imports

import argparse
import Bio
from Bio import SeqIO
import gzip
import pickle
import os
import random


###=============================================================================
### Helper functions

## Component Functions ================

def encode_kmer(kmer):
    char_map = {'A':0, 'C':1, 'G':2, 'T':3}
    kmer_num = 0
    for c in kmer:
        kmer_num = (kmer_num << 2) | char_map[c]
    return kmer_num

def decode_kmer(kmer_num, k=32):
    char_map = {0:'A', 1:'C', 2:'G', 3:'T'}
    kmer = []
    for i in range(k):
        nucleotide_code = kmer_num & 3 # 0x11
        kmer.append(char_map[nucleotide_code])
        kmer_num >>= 2
    return ''.join(reversed(kmer))

def n_check(seq, sample_start, sample_end, skipped_regions):
    "Skip over regions with N"
    gap_index = seq[sample_start:sample_end].find('N')        
    if gap_index > -1:
        sample_end = sample_start + gap_index + 1
        skipped_regions.append((sample_start, sample_end))
        sample_start = sample_end
    return(skipped_regions, sample_start, sample_end)

def sample_scan(seq, start, end, k, seen_kmers, sample_seen_kmers, masked_starts, dedup_retain):
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
            if rng.random() > dedup_retain:
                continue
            else:
                # Handle finding a repeat
                return sample_seen_kmers, masked_starts, True, kmer_start_idx
    return sample_seen_kmers, masked_starts, False, end

## I/O Functions ================

def type_check(file):
    """
    Check if file is gzipped or not
    """
    if file.endswith(".gz") or file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fna") or file.endswith(".txt") or file.endswith(".list"):
        pass
    else:
        raise ValueError("Error: could not determine file type. Supported types are .fasta, .fa, .fna, .fasta.gz, .txt, .list")
        
# Function to check a sample for duplicates, allowing some rate of duplicates
def check_sample(seq, seen_kmers, k, allowed_duplicate_rate):

    # Data structure for new kmers in this sample - see note above
    sample_seen_kmers = set()
    
    # Loop through all kmers in this possible sample
    for kmer_start_idx in range(len(seq)-k+1):

        # Get kmer and its encoding for storage
        kmer = seq[kmer_start_idx:kmer_start_idx+k]
        kmer_num = encode_kmer(kmer)

        # If we haven't seen this kmer before, record it in the sample kmers
        if kmer_num not in seen_kmers and kmer_num not in sample_seen_kmers:
            sample_seen_kmers.add(kmer_num)

        # If we've seen this kmer before, stop analyzing this sample
        else:

            # Stochasticity is built in here here; if random() < allowed duplicate rate, let through this duplicate
            if allowed_duplicate_rate > 0 and rng.random() < allowed_duplicate_rate:
                continue

            # Handle finding a repeat
            sample_end = kmer_start_idx
            return sample_seen_kmers, kmer_start_idx

    return sample_seen_kmers, -1

## Output Functions ================

def writeout_bed(name, regions, bedfile, tpl_input=True, seq_len=None):
    """
    Create bed files for regions or just a list of start sites for beginning positions
    """
    #with open(bedfile, 'w') as f:
    if len(regions) < 1:
        return
    if tpl_input:
        for start, end in regions:
            bedfile.write(f"{name}\t{start}\t{end}\n")
    else:
        for start in regions:
            if seq_len is None:
                bedfile.write(f"{name}\t{start}\n")
            else:
                bedfile.write(f"{name}\t{start}\t{start+seq_len}\n")

def writeout_kmers(seen_kmer_set, file_basename):
    """
    Dump seen kmers into a pickle file that can be used as input for next dedup file
    """
    outfile = file_basename + ".kmers.pkl"
    with open(outfile, 'wb') as f:
        pickle.dump(seen_kmer_set, f)

def output_dump(local_dict, file_basename, k):
    """
    Output all bed files
    """
    with open(file_basename + ".samples.bed", 'w') as sample_regions:
        with open(file_basename + ".masks.bed", 'w') as masked_starts:
            with open(file_basename + ".ignored.bed", 'w') as skipped_regions:
                for seqname in local_dict.keys():
                    writeout_bed(seqname, local_dict[seqname]["sample_regions"], sample_regions)
                    writeout_bed(seqname, local_dict[seqname]["masked_starts"], masked_starts, tpl_input=False, seq_len=k)
                    writeout_bed(seqname, local_dict[seqname]["skipped_regions"], skipped_regions)

## Testing Functions ================

def test_with_fasta(fasta, k, sample_len, seen_kmers):
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

    print(f"Done with {fasta}")


def test_with_sequence(sequence, k, sample_len, seen_kmers):
    sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate(sequence, k, sample_len, seen_kmers=seen_kmers)
    print(sample_regions)
    print(masked_starts)
    print(skipped_regions)
    print(seen_kmers)


###=============================================================================
### Main functions

## Main Deduplication ================
def deduplicate_seq(seq, seen_kmers, args):

    ## Housekeeping ========

    # Collect needed args
    k = args.kmer
    sample_len = args.sample_len
    min_sample_len = args.min_sample_len
    no_overlap_samples = args.no_overlap
    allowed_duplicate_rate = args.retain
    seed = args.seed

    # Instantiate seen_kmers if None
    if seen_kmers is None:
        seen_kmers = set()

    # Set up storage for samples, masked regions, and skipped regions
    sample_regions = []
    masked_starts = []
    skipped_regions = []

    # Check validity of inputs
    if min_sample_len is not None and len(seq) < min_sample_len:
        print("Warning: the sequence length is less than min_sample_len. Skipping this sequence.")
        return sample_regions, masked_starts, skipped_regions, seen_kmers
    if min_sample_len is None and len(seq) < sample_len:
        print("Warning: the sequence length is less than sample_len. Skipping this sequence.")
        return sample_regions, masked_starts, skipped_regions, seen_kmers

    ## Deduplication ========

    # Needed global vars
    max_start_idx = len(seq) - sample_len # final possible start index
    sample_start = 0 # start index of the current sample

    # Investigate every possible sample
    while sample_start <= max_start_idx:

        # Get boundary for this possible sample
        sample_end = sample_start + sample_len

        # Check for Ns
        N_index = seq[sample_start:sample_end].find('N')
        # If we find an N
        if N_index > -1:
            # If the first N occurs before the minimum sample length, there's no valid sample here; skip ahead
            if N_index < min_sample_len:
                # Record entire region before the N as skipped
                if N_index > 0:
                    skipped_regions.append((sample_start, sample_start+N_index))
                sample_start = sample_start + N_index + 1
                skipped_ahead = True
                continue
            # If the first N occurs after the minimum sample length, set the N to be the sample end and continue on to check for the validity of this sample
            else:
                sample_end = sample_start + N_index


        # Check this sample
        sample_seen_kmers, first_repeat_idx = check_sample(seq[sample_start:sample_end], seen_kmers, k, allowed_duplicate_rate)

        # If we didn't find any duplicates:
        # 1. Save the coordinates of this sample as a valid sample
        # 2. Add all kmers from this sample to the global set
        # 3. Set the next sample start based on whether no_overlap is true
        if first_repeat_idx < 0:
            sample_regions.append((sample_start, sample_end))
            seen_kmers.update(sample_seen_kmers)
            sample_start = sample_end if no_overlap_samples else sample_end - k + 1

        # If we did find a duplicate:
        # 1. Add the region from sample start to the duplicate start to the skipped region list
        # 2. Record the duplicate start idx in masked regions
        # 3. Set the next sample start to be 1 + the start position of the duplicate
        else:
            adj_first_repeat_idx = sample_start + first_repeat_idx
            if adj_first_repeat_idx > sample_start:
                skipped_regions.append((sample_start, adj_first_repeat_idx))
            masked_starts.append(adj_first_repeat_idx)
            sample_start = adj_first_repeat_idx + 1

    # Handle possible skipped sequence at the end
    if sample_start < len(seq):
        skipped_regions.append((sample_start, len(seq)))

    return sample_regions, masked_starts, skipped_regions, seen_kmers


def deduplicate_genome(fasta, seen_kmers, args):

    print("Here in deduplicate_genome()")

    # Keep dict associating the regions with the sequence name
    local_dict = {}

    # Set output prefix
    n_suffixes = 2 if fasta.endswith(".gz") else 1
    fasta_basename = '.'.join(os.path.basename(fasta).split('.')[:-(n_suffixes)])
    print(f"fasta basename: {fasta_basename}")
    out_prefix = os.path.join(args.output_dir, fasta_basename)

    # Flexible function to open either a regular or gzipped file
    def open_maybe_gzip(fname):
        with open(fname, "rb") as raw:
            signature = raw.read(2)
        if signature == b"\x1f\x8b":
            return gzip.open(fname, "rt")
        return open(fname, "rt")

    # Read in fasta
    with open_maybe_gzip(fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq)
            seqname = record.id

            print(f"{seqname}: {len(sequence)} bp")

            # Deduplicate this sequence
            sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate_seq(sequence, seen_kmers, args)

            print(f"Sample regions: {sample_regions}")
            print(f"Masked starts: {masked_starts}")
            print(f"Skipped regions: {skipped_regions}")

            # Associate deduplication data with the sequence name
            local_dict[seqname] = {
                "sample_regions": sample_regions,
                "masked_starts": masked_starts,
                "skipped_regions": skipped_regions,
            }

        output_dump(local_dict, out_prefix, args.kmer)

    writeout_kmers(seen_kmers, out_prefix)

    print(f"Done with {fasta}")
    return seen_kmers


def deduplicate(args):

    print("Here in deduplicate()")
    print(f"args: {args}")

    # Read input file of genome locations
    with open(args.input, 'r') as f:
        fastas = [line.rstrip('\n') for line in f.readlines()]

    # Filter down to only the fastas that we can find and issue warning about
    # any we can't find
    valid_fastas = [fasta for fasta in fastas if os.path.isfile(fasta)]
    invalid_fastas = list(set(fastas).difference(set(valid_fastas)))
    if len(invalid_fastas) > 0:
        print("Warning: could not find the following fastas:")
        print(invalid_fastas)

    # Initialize seen kmers here
    if args.seen_kmers is None:
        seen_kmers = set()
    else:
        print("Loading seen kmers")
        seen_kmers = pickle.load(open(args.seen_kmers, "rb"))

    # Iterate over fastas and deduplicate each one
    for fasta in valid_fastas:

        # Deduplicate this fasta
        print(f"Deduplicating {fasta}")
        seen_kmers = deduplicate_genome(fasta, seen_kmers, args)


###=============================================================================
### Call to main

def __main__():

    ## Collect input args
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input list of FASTA files")
    parser.add_argument("-k", "--kmer", type=int, default=32, help="Kmer size (default: 32)")
    parser.add_argument("-l", "--sample_len", type=int, default=1000, help="Sample length (default: 1000)")
    parser.add_argument("-m", "--min_sample_len", type=int, default=None, help="Minimum sample length (default: 50)")
    parser.add_argument("-o", "--output_dir", default="dedup_out", help="Output directory (default: dedup_out/)")
    parser.add_argument("-p", "--seen_kmers", default=None, help="Pickle file containing seen kmers (default: None)")
    parser.add_argument("-r", "--retain", type=float, default=0.0, help="Likelihood a duplicate kmer will be allowed through as a value from [0,1]")
    parser.add_argument("--no-overlap", action="store_true", help="Keep neighboring samples discrete rather than overlapping by k-1 bases")
    parser.add_argument("-seed", "--seed", type=int, default=123, help="Random seed for reproducibility")
    # Hidden test function -- pass in a file (.fa or .txt) here with expected results to verify correctness
    parser.add_argument("-T", "--test", help=argparse.SUPPRESS)
    # Hidden test kmer set input for testing. If not included, will start with empty kmer set
    parser.add_argument("-I", "--test-input-kmers",  help=argparse.SUPPRESS)
    parser.add_argument("-O", "--test-output-kmers", help=argparse.SUPPRESS)
    args = parser.parse_args()

    ## Process input args, checking for validity

    # Input fasta must be a valid file
    if not os.path.isfile(args.input):
        raise("Error: could not find supplied list of fasta files")

    # Kmer must be a positive number
    if args.kmer < 1:
        raise("Error: kmer size must be a positive integer")

    # Set min kmer length to sample length if None
    if args.min_sample_len is None:
        args.min_sample_len = args.sample_len

    # Seen kmers must be either none or a valid pickle file
    if args.seen_kmers is not None:
        if not os.path.isfile(args.seen_kmers):
            raise("Error: could not find supplied seen kmers file")

    # Retain rate must be between 0 and 1
    if args.retain is not None and (args.retain > 1.0 or args.retain < 0.0):
        raise("Error: Retention rate of duplicate kmers cannot be greater than 1.0 or less than 0.0")

    # If sample len is smaller than k, issue warning that no deduplication will occur. Not an error bc this can be a way to create the raw sample sets
    if args.sample_len < args.kmer:
        print("Warning: sample len is smaller than k, meaning no deduplication will occur")

    # If min sample len is higher than sample len, set min to equal the default sample len
    if args.min_sample_len > args.sample_len:
        print("Warning: min sample length cannot be bigger than the default sample length; defaulting to equal the standard sample length")
        args.min_sample_len = args.sample_len

    # If kmer is very small, issue a warning
    if args.kmer < 16:
        print("Warning: small kmer sizes will result in very strict deduplication. Consider increasing the kmer size")

    # Create output directory if it doesn't already exist
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    # Test all input files are valid extension types
    for n in args.input:
        type_check(n)

    
    ## Run deduplication
    if not args.test:
        deduplicate(args)
    else:
        if args.input.endswith(".txt"):
            with open(args.input, 'r') as f:
                sequence = f.read().strip()
            test_with_sequence(sequence, args.kmer, args.sample_len, None)
        else:
            file_basename = args.input().split('.')[0]
            file_outputname=file_basename + ".pickle"
            input_kmers = pickle.load(open(args.test_input_kmers, "rb")) if args.test_input_kmers is not None else set()
            assert args.ouptput_kmers is not None, "Error: must provide output kmer file when testing"
            output_compare = pickle.load(open(args.output_kmers, "rb"))
            test_with_fasta(args.input, args.kmer, args.sample_len, input_kmers)
            assert output_kmers == pickle.load(file_outputname), "Error: output kmers do not match expected output"
            assert

if __name__ == "__main__":
    __main__()

