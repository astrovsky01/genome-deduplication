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
from Bio import SeqIO
import gzip
import json
import pickle
import os
import random as rng
import re
import sys
import tracemalloc


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

def get_fasta_basename(fasta):
    n_suffixes = 2 if fasta.endswith(".gz") else 1
    fasta_basename = '.'.join(os.path.basename(fasta).split('.')[:-(n_suffixes)])
    return fasta_basename

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

# Function to check a sample for duplicates, allowing some rate of duplicates
def check_sample(seq, seen_kmers, k, allowed_duplicate_rate, min_sample_len, no_overlap_samples, N_policy="allow_none"):

    ## Initialize all data structures returned by this function

    # Data structure for new kmers in this sample
    sample_seen_kmers = {}
    
    # Sequence attributes
    sample_end_coord = len(seq)
    duplicate_start_idx = -1
    ambiguous_idx = -1
    ignored_region = None
    next_start_offset = len(seq) if no_overlap_samples else len(seq) - k + 1
    
    # # Check for Ns
    if N_policy == "allow_none":
        N_index = seq.find('N')
        # If we find an N
        if N_index > -1:
            # If the first N occurs before the minimum sample length, there's no valid sample here; handle like a failed sample
            if N_index < min_sample_len:
                # Record position of the N in the list of ambiguous characters
                ambiguous_idx = N_index
                # Record entire region before the N as skipped
                if N_index > 0:
                    #NOT ADDING K HERE
                    ignored_region = (0, N_index)
                next_start_offset = N_index + 1
                sample_end_coord = -1
            else:
                # Set N index to be the end of this possible sample
                sample_end_coord = N_index

    # Loop through all kmers in this possible sample
    if sample_end_coord >= min_sample_len:
        for kmer_start_idx in range(sample_end_coord-k+1):

            # Get kmer and its encoding for storage
            kmer = seq[kmer_start_idx:kmer_start_idx+k]
            kmer_num = encode_kmer(kmer)

            # If we haven't seen this kmer before, record it in the sample kmers
            if kmer_num not in seen_kmers and kmer_num not in sample_seen_kmers:
                sample_seen_kmers[kmer_num] = kmer_start_idx

            # If we've seen this kmer before, stop analyzing this sample
            else:

                # Opportunity to process this duplicate with some function to determine whether to ignore it
                # Simple chance, function based on previous occurrence count, etc.
                if allowed_duplicate_rate > 0 and rng.random() < allowed_duplicate_rate:
                    continue

                # At this point, we have found a definite duplicate

                # If internal repeat (offender is a repeat of original):
                    # If offender comes before min_sample_len:
                        # Although we don't save sample seen kmers, we know that starting again before the original kmer will run into the same repeat
                        # Therefore, next sample needs to start after the original
                        # -> Invalid sample
                        # -> Ignored region from start to original + 1
                        # -> No duplicate found
                        # -> Next start = original + 1
                    # else:
                        # We will save this sample and add its kmers to the global set, so any new sample before the offender will be invalidated by the original
                        # Therefore, next sample needs to start after the offender
                        # -> valid sample up to offender position
                        # -> Next start = offender + 1
                # else (repeat matches the global, not the local, set)
                    # If offender comes before min_sample_len:
                        # -> Invalid sample
                        # -> Next start = offender + 1
                    # else:
                        # -> valid sample up to offender position
                        # -> Next start = offender + 1

                # If the offending kmer is within min_sample_len
                if kmer_start_idx < min_sample_len:

                    # This is an invalid sample
                    sample_end_coord = -1

                    # Check if this is an internal duplicate
                    repeat_match_idx = sample_seen_kmers.get(kmer_num, -1)

                    # If an internal match
                    if repeat_match_idx >= 0:

                        # The next sample should start at the original kmer's position + 1
                        next_start_offset = repeat_match_idx + 1

                        # Record the entire region up to and including the original kmer as ignored
                        # Importantly, we do not denote the original kmer as a duplicate, since it is both not in the global set and is not added to a sample
                        #if repeat_match_idx > 0:
                        ##ADDED K HERE
                        ignored_region = (0, repeat_match_idx+k)

                    # Else if a match to a global kmer
                    else:

                        # The next sample should start at the current position + 1
                        next_start_offset = kmer_start_idx + 1

                        # Record the current kmer as a duplicate
                        duplicate_start_idx = kmer_start_idx

                        # Record the entire region up to but not including the offending kmer as ignored
                        # Include full final kmer position
                        if kmer_start_idx > 0:
                            ##ADDED K HERE
                            ignored_region = (0, kmer_start_idx+k-1)

                    # Since this sample is invalid, we can discard its seen kmers
                    sample_seen_kmers = None

                # If the offending kmer is past min_sample_len
                # Whether the repeat is internal or global, the current kmer is the offending one
                else:

                    # The current kmer is the end boundary of the sample
                    sample_end_coord = kmer_start_idx

                    # The current kmer is a duplicate
                    duplicate_start_idx = kmer_start_idx

                    # The next sample should start at 1 + the current position
                    next_start_offset = kmer_start_idx + 1
                
                # We have finished processing this sample
                break

    return sample_end_coord, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers

# Group masked indices into contiguous bed regions
# e.g. [2,3,4,7,8,20] -> [(2,5), (7,9), (20,21)]
def condense_masked_regions(masked, k=32):
    masked_regions = []
    if len(masked) == 0:
        return masked_regions
    region_start = masked[0]
    for i in range(len(masked)-1):
        # If the next masked kmer is contiguous with the current one
        if masked[i] + 1 != masked[i+1]:
            masked_regions.append((region_start, masked[i] + k))
            region_start = masked[i+1]
    masked_regions.append((region_start, masked[-1] + k))
    return masked_regions

def condense_ambiguous_regions(ambiguous):
    ambiguous_regions = []
    if len(ambiguous) == 0:
        return ambiguous_regions
    region_start = ambiguous[0]
    for i in range(len(ambiguous)-1):
        if ambiguous[i] + 1 != ambiguous[i+1]:
            ambiguous_regions.append((region_start, ambiguous[i] + 1))
            region_start = ambiguous[i+1]
    ambiguous_regions.append((region_start, ambiguous[-1] + 1))
    return ambiguous_regions

## I/O Functions ================

def type_check(file):
    """
    Check if file is gzipped or not
    """
    if file.endswith(".gz") or file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fna") or file.endswith(".txt") or file.endswith(".list"):
        pass
    else:
        raise ValueError("Error: could not determine file type. Supported types are .fasta, .fa, .fna, .fasta.gz, .txt, .list")

# Flexible function to open either a regular or gzipped file
def open_maybe_gzip(fname):
    with open(fname, "rb") as raw:
        signature = raw.read(2)
    if signature == b"\x1f\x8b":
        return gzip.open(fname, "rt")
    return open(fname, "rt")

# Extra term is a value that will be written to every row of the 4th column
def writeout_bed(name, regions, bedfile, extra_col_val=None):
    """
    Create bed files for regions
    """
    if len(regions) < 1:
        return
    extra_term = "" if extra_col_val is None else f"\t{extra_col_val}"
    for start,end in regions:
        bedfile.write(f"{name}\t{start}\t{end}{extra_term}\n")

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
        with open(file_basename + ".masks.bed", 'w') as masked_regions:
            with open(file_basename + ".ignored.bed", 'w') as skipped_regions:
                with open(file_basename + ".ambiguous.bed", 'w') as ambiguous_regions:
                    for seqname in local_dict.keys():
                        writeout_bed(seqname, local_dict[seqname]["sample_regions"], sample_regions)
                        writeout_bed(seqname, local_dict[seqname]["masked_regions"], masked_regions)
                        writeout_bed(seqname, local_dict[seqname]["skipped_regions"], skipped_regions)
                        writeout_bed(seqname, local_dict[seqname]["ambiguous_regions"], ambiguous_regions)

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


def test_with_sequence(args):
    with open(args.input[0], 'r') as f:
        seq = f.read().strip()
    sample_regions, masked_starts, skipped_regions, seen_kmers = deduplicate_seq(seq, args.seen_kmers, args)
    print("Sample regions: ", sample_regions)
    print("Masked starts: ", masked_starts)
    print("Skipped regions: ", skipped_regions)
    print("Seen kmers: ", seen_kmers)
    return sample_regions, masked_starts, skipped_regions, seen_kmers


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
    ambiguous_positions = []

    # Check validity of inputs
    if min_sample_len is not None and len(seq) < min_sample_len:
        print("Warning: the sequence length is less than min_sample_len. Skipping this sequence.")
        return [], [], [(0,len(seq))], [], seen_kmers 
    if min_sample_len is None and len(seq) < sample_len:
        print("Warning: the sequence length is less than sample_len. Skipping this sequence.")
        return [], [], [(0,len(seq))], [], seen_kmers 

    ## Deduplication ========

    # Needed global vars
    max_start_idx = len(seq) - min_sample_len # final possible start index
    sample_start = 0 # start index of the current sample

    # Investigate every possible sample
    while sample_start <= max_start_idx:

        # Get boundary for this possible sample
        sample_end = min(len(seq), sample_start + sample_len) # Checking against seq len is only necessary for final sample

        # Check this sample
        #checked_sample_len, duplicate_start_idx, next_start_offset, ignored_region, sample_seen_kmers = \
        checked_sample_len, duplicate_start_idx, ambiguous_idx, ignored_region, next_start_offset, sample_seen_kmers = \
                check_sample(seq[sample_start:sample_end], seen_kmers, k, allowed_duplicate_rate, min_sample_len, no_overlap_samples)

        #if verbose:
        #    print(f"{sample_start}: checked_sample_len={checked_sample_len}, duplicate_start_idx={duplicate_start_idx}, next_start_offset={next_start_offset}, ignored_region={ignored_region}")

        # Add sample to samples if valid
        if checked_sample_len > -1:
            sample_regions.append((sample_start, sample_start + checked_sample_len))
        
        # Record the duplicate if any was found
        if duplicate_start_idx > -1:
            masked_starts.append(sample_start + duplicate_start_idx)

        # Record the ambiguous index if any was found
        if ambiguous_idx > -1:
            ambiguous_positions.append(sample_start + ambiguous_idx)

        # Add ignored region if any was found
        if ignored_region is not None:
            skipped_regions.append((sample_start+ignored_region[0], sample_start+ignored_region[1]))

        # Update seen kmers with any sample-specific kmers
        if sample_seen_kmers is not None:
            seen_kmers.update(sample_seen_kmers)

        # Set start index of next iteration of the loop
        sample_start = sample_start + next_start_offset

    # Handle possible skipped sequence at the end
    if sample_start < len(seq):
        ## N Check final region and divide between skipped and ambiguous
        final_sequence=seq[sample_start:len(seq)]
        while len(final_sequence) > 0:
            N_index = final_sequence.find('N')
            # If we find an N
            if N_index > -1:
                if N_index != 0:
                    skipped_regions.append((sample_start, sample_start+N_index))
                ambiguous_positions.append(sample_start + N_index)
                sample_start = sample_start + N_index + 1
            else:
                skipped_regions.append((sample_start, len(seq)))
            final_sequence=final_sequence[N_index+1:] if N_index > -1 else ""

    # Convert masked starting indices and ambiguous positions to regions for more condensed bed files
    masked_regions = condense_masked_regions(masked_starts)
    ambiguous_regions = condense_ambiguous_regions(ambiguous_positions)

    return sample_regions, masked_regions, skipped_regions, ambiguous_regions, seen_kmers


def deduplicate_genome(fasta, seen_kmers, save_kmers_to_file, args):

    #print("Here in deduplicate_genome()")

    # Keep dict associating the regions with the sequence name
    local_dict = {}

    # Set output prefix
    fasta_basename = get_fasta_basename(fasta)
    out_prefix = os.path.join(args.output_dir, fasta_basename)

    # Read in fasta
    with open_maybe_gzip(fasta) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequence = str(record.seq)
            seqname = record.id

            # Convert all chars to uppercase and replace ambiguous chars with N
            ambiguous_chars_regex = re.compile(r'[^ACGTN]')
            clean_sequence = ambiguous_chars_regex.sub('N', str(sequence.upper()))

            #print(f"{seqname}: {len(sequence)} bp")

            # Deduplicate this sequence
            sample_regions, masked_regions, skipped_regions, ambiguous_regions, seen_kmers = deduplicate_seq(clean_sequence, seen_kmers, args)

            #print(f"Sample regions: {sample_regions}")
            #print(f"Masked regions: {masked_regions}")
            #print(f"Skipped regions: {skipped_regions}")

            # Associate deduplication data with the sequence name
            local_dict[seqname] = {
                "sample_regions": sample_regions,
                "masked_regions": masked_regions,
                "skipped_regions": skipped_regions,
                "ambiguous_regions": ambiguous_regions
            }

        output_dump(local_dict, out_prefix, args.kmer)

    if save_kmers_to_file:
        writeout_kmers(seen_kmers, out_prefix)

    #print(f"Done with {fasta}")
    return seen_kmers


def deduplicate(args):

    #print("Here in deduplicate()")
    print(f"args: {args}")
    if not args.yes:
    # Create output directory if it doesn't already exist
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)
        else:
            delete_check=input("Output directory already exists. Continue and potentially overwrite files? (y/n): ")
            if delete_check.lower() != 'y' and delete_check.lower() != 'yes':
                print("Exiting...")
                sys.exit(1)
    else:
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)

    # Write (processed) input args to file for reproducibility
    with open(os.path.join(args.output_dir, "config.json"), 'w') as f:
        json.dump(vars(args), f, indent=4)

    # Read input file of genome locations
    # Single fasta input or series of individual fasta entries
    if len(args.input) > 1 or args.input[0].endswith(('.fa', '.fasta', '.fasta.gz', '.fna', '.fna.gz')):
        fastas = args.input
    # File with a list of fasta files as input
    elif len(args.input) == 1 and args.input[0].endswith(('.txt', '.list')):
        with open(args.input[0], 'r') as f:
            fastas = [line.rstrip('\n') for line in f.readlines()]

    # Filter down to only the fastas that we can find and issue warning about
    # any we can't find
    valid_fastas = [fasta for fasta in fastas if os.path.isfile(fasta)]
    invalid_fastas = list(set(fastas).difference(set(valid_fastas)))
    if len(invalid_fastas) > 0:
        print("Warning: could not find the following fastas:")
        print(invalid_fastas)

    # Write basename to file map for all valid fastas
    basename_fasta_file = os.path.join(args.output_dir, "basename_fasta_match.txt")
    with open(basename_fasta_file, 'w') as f:
        for fasta in valid_fastas:
            fasta_basename = get_fasta_basename(fasta)
            f.write(f"{fasta_basename}\t{fasta}\n")

    # Initialize seen kmers here
    if args.seen_kmers is None:
        seen_kmers = set()
    else:
        print("Loading seen kmers")
        seen_kmers = pickle.load(open(args.seen_kmers, "rb"))

    # Iterate over fastas and deduplicate each one
    for i,fasta in enumerate(valid_fastas):

        # Deduplicate this fasta
        print(f"Deduplicating {fasta}")
        save_kmers_to_file = args.save_every > 0 and (i+1) % args.save_every == 0
        seen_kmers = deduplicate_genome(fasta, seen_kmers, save_kmers_to_file, args)

    # Optionally save seen kmers at the end
    if not args.no_save_kmers_at_end:
        final_seen_kmer_file_basename = "final"
        idx_suffix = 1
        while os.path.exists(os.path.join(args.output_dir, f"{final_seen_kmer_file_basename}.kmers.pkl")):
            final_seen_kmer_file_basename = f"final_{idx_suffix}"
            idx_suffix += 1
        writeout_kmers(seen_kmers, os.path.join(args.output_dir, final_seen_kmer_file_basename))


###=============================================================================
### Call to main

def __main__():

    tracemalloc.start()

    ## Collect input args
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="+", help="Input list of FASTA files or a txt file with one FASTA file per line")
    parser.add_argument("-k", "--kmer", type=int, default=32, help="Kmer size (default: 32)")
    parser.add_argument("-l", "--sample_len", type=int, default=1000, help="Sample length (default: 1000)")
    parser.add_argument("-m", "--min_sample_len", type=int, default=None, help="Minimum sample length (default: 50)")
    parser.add_argument("-o", "--output_dir", default="dedup_out", help="Output directory (default: dedup_out/)")
    parser.add_argument("-p", "--seen_kmers", default=None, help="Pickle file containing seen kmers (default: None)")
    parser.add_argument("-r", "--retain", type=float, default=0.0, help="Likelihood a duplicate kmer will be allowed through as a value from [0,1]")
    parser.add_argument("-s", "--save_every", type=int, default=0, help="Save seen kmer set every n samples (default: 0, don't save any before the end)")
    parser.add_argument("--no-overlap", action="store_true", help="Keep neighboring samples discrete rather than overlapping by k-1 bases")
    parser.add_argument("--no-save-kmers-at-end", action="store_true", help="Opt out of saving the seen kmer set at the end of program execution")
    parser.add_argument("-seed", "--seed", type=int, default=123, help="Random seed for reproducibility")
    parser.add_argument("-y", "--yes", action="store_true", help="Automatic overwrite of output directory if it already exists")
    # Hidden test function -- pass in a file (.fa or .txt) here with expected results to verify correctness
    parser.add_argument("-T", "--test", action="store_true", help=argparse.SUPPRESS)
    # Hidden test kmer set input for testing. If not included, will start with empty kmer set
    parser.add_argument("-I", "--test-input-kmers", type=str, help=argparse.SUPPRESS)
    parser.add_argument("-O", "--test-output-kmers", type=str, help=argparse.SUPPRESS)
    args = parser.parse_args()

    ## Process input args, checking for validity

    # Input fasta must be a valid file
    for f in args.input:
        if not os.path.isfile(f):
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

    ## Run deduplication
    print(args)
    if not args.test:
        print("Running deduplication")
        deduplicate(args)
    else:
        input=args.input[0]
        if input.endswith(".txt"):
            with open(input) as f:
                if f.readline().endswith(".fa\n") or f.readline().endswith(".fasta\n") or f.readline().endswith(".fna\n") or f.readline().endswith(".fa.gz\n") or f.readline().endswith(".fasta.gz\n") or f.readline().endswith(".fna.gz\n"):
                    print("Testing with fasta")
                    test_with_fasta(args)
            print("Testing with sequence")
            test_with_sequence(args)
        else:
            print("Testing with fasta")
            file_basename = input().split('.')[0]
            file_outputname=file_basename + ".pickle"
            input_kmers = pickle.load(open(args.test_input_kmers, "rb")) if args.test_input_kmers is not None else set()
            assert args.ouptput_kmers is not None, "Error: must provide output kmer file when testing"
            output_compare = pickle.load(open(args.output_kmers, "rb"))
            test_with_fasta(input, args.kmer, args.sample_len, input_kmers)
            assert output_compare == pickle.load(file_outputname), "Error: output kmers do not match expected output"

if __name__ == "__main__":
    __main__()

