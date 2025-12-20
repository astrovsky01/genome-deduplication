import argparse

parser = argparse.ArgumentParser(description="Identify unused kmers from deduplicated kmer list.")
parser.add_argument("-d", "--deduped_kmers", required=True, help="Path to deduplicated kmers file.")
parser.add_argument("-f", "--ignored_file", required=True, help="Ignored bed file containing kmers to check")
parser.add_argument("-n", "--file_name", required=True, help="Source file basename")
parser.add_argument("-k", "--kmer_size", type=int, required=True, help="Size of kmers.")
parser.add_argument("-o", "--outfile", required=True, help="Output file for unused kmers.")
args = parser.parse_args()

all_deduped_kmers = args.deduped_kmers
ignored_file = args.ignored_file
kmer_size = args.kmer_size
file_name = args.file_name
output_file = args.outfile

# Load deduped kmers into hash table for fast lookup
kmer_set = set()
counter = 0
collisions = 0
print("Loading deduped kmers into hash table...")
with open(all_deduped_kmers, 'r') as f:
    for line in f:
        kmer = line.strip().split()[0].upper()
        if kmer in kmer_set:
            collisions += 1
        kmer_set.add(kmer)
        counter += 1
        if counter % 200000 == 0:
            print(f"Loaded {counter} kmers...")
print(f"Loaded {len(kmer_set)} unique kmers from {counter} total lines with {collisions} collisions")


with open(ignored_file, 'r') as f, open(output_file, 'a') as out_f:
    for line in f:
        # breakdown bed file line to seqname, start, end, sequence
        name_coords, sequence = line.strip().split('\t')
        name = name_coords.split(':')[0]
        coords = name_coords.split(':')[1]
        start = int(coords.split('-')[0])
        end = int(coords.split('-')[1])
        if len(sequence) < kmer_size:
            continue
            
        # Find consecutive kmers that were never used to reduce file size. Merge kmers into longer sequences
        current_seq = ""
        current_start = None
        current_end = None
        
        for i in range(len(sequence) - kmer_size + 1):
            kmer = sequence[i:i + kmer_size].upper()
            kmer_start = start + i
            kmer_end = start + i + kmer_size - 1
            
            if kmer not in kmer_set:
                if current_seq == "":
                    # Start new sequence
                    current_seq = kmer
                    current_start = kmer_start
                    current_end = kmer_end
                else:
                    # Extend current sequence (next consecutive k-mer)
                    current_seq += kmer[-1]  # Add only the last nucleotide
                    current_end = kmer_end
            else:
                print("FOUND ONE")
                print(kmer, "at", kmer_start, "to", kmer_end, "in", name)
                # K-mer exists in set - write any accumulated sequence
                if current_seq != "":
                    out_f.write(f"{name}\t{current_start}\t{current_end}\t{current_seq}\t{file_name}\n")
                current_seq = ""
                current_start = None
                current_end = None
        
        # Write final sequence if any
        if current_seq != "":
            out_f.write(f"{name}\t{current_start}\t{current_end}\t{current_seq}\t{file_name}\n")
        

