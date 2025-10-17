/**
 * Genome Deduplication Tool - C++ Implementation
 * 
 * Deduplicates and samples genomic sequences simultaneously
 * Only considers a kmer if it is within a sample
 * 
 * Basic algorithm:
 *     1. Sample sequence to first `sample_len` bases
 *     2. Iterate over all kmers in sample
 *     3. If kmer not in seen_kmers, add it to seen_kmers
 *     4. If kmer in seen_kmers, jump to next possible sample and repeat
 * 
 * The program divides a fasta into three categories:
 *     1. Samples - valid regions with no repeats, on which we train
 *     2. Masked regions - regions containing duplicates
 *     3. Ignored regions - regions that we threw out because we found a
 *        duplicate later on in the sample region
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <random>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <cstring>
#include <cctype>
#include <zlib.h>

namespace fs = std::filesystem;

// Structure to hold genomic regions (start inclusive, end exclusive)
struct Region {
    size_t start;
    size_t end;
    
    Region(size_t s, size_t e) : start(s), end(e) {}
};

// Structure to hold deduplication results for a sequence
struct DeduplicationResult {
    std::vector<Region> sample_regions;
    std::vector<Region> masked_regions;
    std::vector<Region> skipped_regions;
    std::vector<Region> ambiguous_regions;
};

// Structure to hold program arguments
struct Args {
    std::vector<std::string> input;
    int kmer = 32;
    int sample_len = 1000;
    int min_sample_len = -1; // Will be set to sample_len if not specified
    std::string output_dir = "dedup_out";
    std::string seen_kmers_file;
    double retain = 0.0;
    int save_every = 0;
    bool no_overlap = false;
    bool no_save_kmers_at_end = false;
    int seed = 123;
};

// Global random number generator
std::mt19937 rng_gen;

/**
 * Encode a k-mer string into a 64-bit integer
 * A=0, C=1, G=2, T=3
 */
uint64_t encode_kmer(const std::string& kmer) {
    uint64_t kmer_num = 0;
    for (char c : kmer) {
        kmer_num = (kmer_num << 2);
        switch(c) {
            case 'A': kmer_num |= 0; break;
            case 'C': kmer_num |= 1; break;
            case 'G': kmer_num |= 2; break;
            case 'T': kmer_num |= 3; break;
            default: 
                throw std::runtime_error("Invalid nucleotide in k-mer: " + std::string(1, c));
        }
    }
    return kmer_num;
}

/**
 * Decode a k-mer integer back to string
 */
std::string decode_kmer(uint64_t kmer_num, int k = 32) {
    std::string kmer(k, 'A');
    for (int i = k - 1; i >= 0; --i) {
        int nucleotide_code = kmer_num & 3;
        switch(nucleotide_code) {
            case 0: kmer[i] = 'A'; break;
            case 1: kmer[i] = 'C'; break;
            case 2: kmer[i] = 'G'; break;
            case 3: kmer[i] = 'T'; break;
        }
        kmer_num >>= 2;
    }
    return kmer;
}

/**
 * Get the basename of a fasta file (without path and extensions)
 */
std::string get_fasta_basename(const std::string& fasta) {
    fs::path p(fasta);
    std::string filename = p.filename().string();
    
    // Remove .gz if present
    if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
        filename = filename.substr(0, filename.size() - 3);
    }
    
    // Remove fasta extensions
    size_t last_dot = filename.find_last_of('.');
    if (last_dot != std::string::npos) {
        return filename.substr(0, last_dot);
    }
    return filename;
}

/**
 * Condense a list of indices into contiguous regions
 * e.g. [2,3,4,7,8,20] -> [(2,5), (7,9), (20,21)]
 */
std::vector<Region> condense_masked_regions(const std::vector<size_t>& masked) {
    std::vector<Region> masked_regions;
    if (masked.empty()) {
        return masked_regions;
    }
    
    size_t region_start = masked[0];
    for (size_t i = 0; i < masked.size() - 1; ++i) {
        if (masked[i] + 1 != masked[i + 1]) {
            masked_regions.emplace_back(region_start, masked[i] + 1);
            region_start = masked[i + 1];
        }
    }
    masked_regions.emplace_back(region_start, masked.back() + 1);
    return masked_regions;
}

/**
 * Check a sample for duplicates
 */
struct CheckSampleResult {
    int checked_sample_len;
    int duplicate_start_idx;
    int ambiguous_idx;
    std::pair<int, int> ignored_region; // -1, -1 if none
    size_t next_start_offset;
    std::unordered_map<uint64_t, size_t> sample_seen_kmers;
    bool valid_sample_kmers;
};

CheckSampleResult check_sample(
    const std::string& seq,
    const std::unordered_set<uint64_t>& seen_kmers,
    int k,
    double allowed_duplicate_rate,
    int min_sample_len,
    bool no_overlap_samples,
    const std::string& N_policy = "allow_none",
    bool debug = false
) {
    CheckSampleResult result;
    result.checked_sample_len = -1;
    result.duplicate_start_idx = -1;
    result.ambiguous_idx = -1;
    result.ignored_region = {-1, -1};
    result.valid_sample_kmers = false;
    
    size_t sample_end_coord = seq.length();
    result.next_start_offset = no_overlap_samples ? seq.length() : seq.length() - k + 1;
    
    // Check for Ns
    if (N_policy == "allow_none") {
        size_t N_index = seq.find('N');
        if (N_index != std::string::npos) {
            if (N_index < static_cast<size_t>(min_sample_len)) {
                result.ambiguous_idx = N_index;
                if (N_index > 0) {
                    result.ignored_region = {0, static_cast<int>(N_index)};
                }
                result.next_start_offset = N_index + 1;
                result.checked_sample_len = -1;
                return result;
            } else {
                sample_end_coord = N_index;
            }
        }
    }
    
    // Loop through all kmers in this possible sample
    if (sample_end_coord >= static_cast<size_t>(min_sample_len)) {
        for (size_t kmer_start_idx = 0; kmer_start_idx <= sample_end_coord - k; ++kmer_start_idx) {
            std::string kmer = seq.substr(kmer_start_idx, k);
            uint64_t kmer_num = encode_kmer(kmer);
            
            if (debug && kmer_start_idx == 0) {
                std::cout << "    First kmer: " << kmer << std::endl;
                std::cout << "    Encoded: " << kmer_num << std::endl;
                std::cout << "    In seen_kmers: " << (seen_kmers.find(kmer_num) != seen_kmers.end()) << std::endl;
            }
            
            // If we haven't seen this kmer before, record it
            if (seen_kmers.find(kmer_num) == seen_kmers.end() &&
                result.sample_seen_kmers.find(kmer_num) == result.sample_seen_kmers.end()) {
                result.sample_seen_kmers[kmer_num] = kmer_start_idx;
            } else {
                if (debug) {
                    std::cout << "    Found dup at kmer_start_idx=" << kmer_start_idx 
                              << " kmer=" << kmer.substr(0,8) << "..."
                              << " in_seen=" << (seen_kmers.find(kmer_num) != seen_kmers.end())
                              << " in_sample=" << (result.sample_seen_kmers.find(kmer_num) != result.sample_seen_kmers.end())
                              << " sample_size=" << result.sample_seen_kmers.size()
                              << std::endl;
                    
                    // Check if internal
                    auto it_internal = result.sample_seen_kmers.find(kmer_num);
                    if (it_internal != result.sample_seen_kmers.end()) {
                        std::cout << "      -> INTERNAL match, original at " << it_internal->second << std::endl;
                    }
                }
                // Found a duplicate - check if we should allow it
                if (allowed_duplicate_rate > 0) {
                    std::uniform_real_distribution<double> dist(0.0, 1.0);
                    if (dist(rng_gen) < allowed_duplicate_rate) {
                        continue;
                    }
                }
                
                // Process this duplicate
                if (kmer_start_idx < static_cast<size_t>(min_sample_len)) {
                    result.checked_sample_len = -1;
                    
                    // Check if internal duplicate
                    auto it = result.sample_seen_kmers.find(kmer_num);
                    if (it != result.sample_seen_kmers.end()) {
                        // Internal match
                        result.next_start_offset = it->second + 1;
                        result.ignored_region = {0, static_cast<int>(it->second + 1)};
                    } else {
                        // Global match
                        result.next_start_offset = kmer_start_idx + 1;
                        result.duplicate_start_idx = kmer_start_idx;
                        if (kmer_start_idx > 0) {
                            result.ignored_region = {0, static_cast<int>(kmer_start_idx)};
                        }
                    }
                    result.valid_sample_kmers = false;
                } else {
                    // Offending kmer past min_sample_len
                    result.checked_sample_len = kmer_start_idx;
                    result.duplicate_start_idx = kmer_start_idx;
                    result.next_start_offset = kmer_start_idx + 1;
                    result.valid_sample_kmers = true;
                }
                break;
            }
        }
        
        // If no duplicate found (valid_sample_kmers is still false from initialization if no duplicate was found)
        if (result.checked_sample_len == -1 && result.duplicate_start_idx == -1 && result.valid_sample_kmers == false && result.ignored_region.first == -1) {
            result.checked_sample_len = sample_end_coord;
            result.valid_sample_kmers = true;
        }
    }
    
    return result;
}

/**
 * Deduplicate a single sequence
 */
DeduplicationResult deduplicate_seq(
    const std::string& seq,
    std::unordered_set<uint64_t>& seen_kmers,
    const Args& args,
    bool debug = false
) {
    DeduplicationResult result;
    
    // Check validity of inputs
    if (seq.length() < static_cast<size_t>(args.min_sample_len)) {
        std::cerr << "Warning: sequence length is less than min_sample_len. Skipping this sequence." << std::endl;
        result.skipped_regions.emplace_back(0, seq.length());
        return result;
    }
    
    std::vector<size_t> masked_starts;
    std::vector<size_t> ambiguous_positions;
    
    size_t max_start_idx = seq.length() - args.min_sample_len;
    size_t sample_start = 0;
    
    // Investigate every possible sample
    while (sample_start <= max_start_idx) {
        size_t sample_end = std::min(seq.length(), sample_start + args.sample_len);
        
        std::string sample_seq = seq.substr(sample_start, sample_end - sample_start);
        
        CheckSampleResult check_result = check_sample(
            sample_seq,
            seen_kmers,
            args.kmer,
            args.retain,
            args.min_sample_len,
            args.no_overlap,
            "allow_none",
            debug && sample_start < 200
        );
        
        if (debug && sample_start < 200) {
            std::cout << "  sample_start=" << sample_start 
                      << " checked_len=" << check_result.checked_sample_len
                      << " dup_idx=" << check_result.duplicate_start_idx
                      << " ignored=(" << check_result.ignored_region.first << "," << check_result.ignored_region.second << ")"
                      << " next_offset=" << check_result.next_start_offset
                      << " sample_kmers=" << check_result.sample_seen_kmers.size()
                      << std::endl;
        }
        
        // Add sample to samples if valid
        if (check_result.checked_sample_len > -1) {
            result.sample_regions.emplace_back(
                sample_start,
                sample_start + check_result.checked_sample_len
            );
        }
        
        // Record the duplicate if any was found
        if (check_result.duplicate_start_idx > -1) {
            masked_starts.push_back(sample_start + check_result.duplicate_start_idx);
        }
        
        // Record the ambiguous index if any was found
        if (check_result.ambiguous_idx > -1) {
            ambiguous_positions.push_back(sample_start + check_result.ambiguous_idx);
        }
        
        // Add ignored region if any was found
        if (check_result.ignored_region.first != -1) {
            result.skipped_regions.emplace_back(
                sample_start + check_result.ignored_region.first,
                sample_start + check_result.ignored_region.second
            );
        }
        
        // Update seen kmers with sample-specific kmers
        if (check_result.valid_sample_kmers) {
            for (const auto& kv : check_result.sample_seen_kmers) {
                seen_kmers.insert(kv.first);
            }
        }
        
        // Set start index of next iteration
        sample_start = sample_start + check_result.next_start_offset;
    }
    
    // Handle possible skipped sequence at the end
    if (sample_start < seq.length()) {
        result.skipped_regions.emplace_back(sample_start, seq.length());
    }
    
    // Convert masked starting indices and ambiguous positions to regions
    result.masked_regions = condense_masked_regions(masked_starts);
    result.ambiguous_regions = condense_masked_regions(ambiguous_positions);
    
    return result;
}

/**
 * Write out BED file regions
 */
void writeout_bed(
    const std::string& name,
    const std::vector<Region>& regions,
    std::ofstream& bedfile
) {
    for (const auto& region : regions) {
        bedfile << name << "\t" << region.start << "\t" << region.end << "\n";
    }
}

/**
 * Output all bed files
 */
void output_dump(
    const std::vector<std::pair<std::string, DeduplicationResult>>& local_dict,
    const std::string& file_basename
) {
    std::ofstream sample_regions(file_basename + ".samples.bed");
    std::ofstream masked_regions(file_basename + ".masks.bed");
    std::ofstream skipped_regions(file_basename + ".ignored.bed");
    std::ofstream ambiguous_regions(file_basename + ".ambiguous.bed");
    
    for (const auto& entry : local_dict) {
        const std::string& seqname = entry.first;
        const DeduplicationResult& result = entry.second;
        
        writeout_bed(seqname, result.sample_regions, sample_regions);
        writeout_bed(seqname, result.masked_regions, masked_regions);
        writeout_bed(seqname, result.skipped_regions, skipped_regions);
        writeout_bed(seqname, result.ambiguous_regions, ambiguous_regions);
    }
}

/**
 * Simple FASTA parser with gzip support
 */
class FastaParser {
private:
    gzFile gz_file;
    std::string current_name;
    std::string current_seq;
    bool has_next;
    bool is_gzipped;
    std::string pending_line;
    bool has_pending;
    
    bool getline(std::string& line) {
        if (has_pending) {
            line = pending_line;
            has_pending = false;
            return true;
        }
        
        line.clear();
        char buffer[4096];
        
        if (gzgets(gz_file, buffer, sizeof(buffer)) != nullptr) {
            line = buffer;
            // Remove trailing newline
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        }
        return false;
    }
    
    void putback_line(const std::string& line) {
        pending_line = line;
        has_pending = true;
    }
    
public:
    FastaParser(const std::string& filename) : has_next(true), has_pending(false) {
        // Check if file is gzipped
        is_gzipped = (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz");
        
        gz_file = gzopen(filename.c_str(), "rb");
        if (gz_file == nullptr) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        advance();
    }
    
    ~FastaParser() {
        if (gz_file != nullptr) {
            gzclose(gz_file);
        }
    }
    
    void advance() {
        current_seq.clear();
        current_name.clear();
        
        std::string line;
        bool found_header = false;
        
        while (getline(line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                if (found_header) {
                    // Put the line back for next iteration
                    putback_line(line);
                    has_next = true;
                    return;
                }
                current_name = line.substr(1);
                // Extract just the ID (first word)
                size_t space_pos = current_name.find(' ');
                if (space_pos != std::string::npos) {
                    current_name = current_name.substr(0, space_pos);
                }
                found_header = true;
            } else if (found_header) {
                current_seq += line;
            }
        }
        
        has_next = found_header && !current_seq.empty();
    }
    
    bool hasNext() const {
        return has_next;
    }
    
    std::string getName() const {
        return current_name;
    }
    
    std::string getSequence() const {
        return current_seq;
    }
};

/**
 * Clean sequence - convert to uppercase and replace ambiguous chars with N
 */
std::string clean_sequence(const std::string& seq) {
    std::string cleaned;
    cleaned.reserve(seq.length());
    
    for (char c : seq) {
        char upper_c = std::toupper(c);
        if (upper_c == 'A' || upper_c == 'C' || upper_c == 'G' || upper_c == 'T' || upper_c == 'N') {
            cleaned += upper_c;
        } else {
            cleaned += 'N';
        }
    }
    
    return cleaned;
}

/**
 * Deduplicate a genome file
 */
void deduplicate_genome(
    const std::string& fasta,
    std::unordered_set<uint64_t>& seen_kmers,
    bool save_kmers_to_file,
    const Args& args
) {
    // Use vector of pairs to maintain insertion order
    std::vector<std::pair<std::string, DeduplicationResult>> local_dict;
    
    std::string fasta_basename = get_fasta_basename(fasta);
    std::string out_prefix = args.output_dir + "/" + fasta_basename;
    
    FastaParser parser(fasta);
    
    while (parser.hasNext()) {
        std::string seqname = parser.getName();
        std::string sequence = parser.getSequence();
        
        std::string clean_seq = clean_sequence(sequence);
        
        //std::cout << seqname << ": " << sequence.length() << " bp" << std::endl;
        
        DeduplicationResult result = deduplicate_seq(clean_seq, seen_kmers, args);
        
        // Match Python's sys.getsizeof - rough estimate of memory usage
        // Python's sys.getsizeof includes overhead, approximate with count * 8 bytes per entry + overhead
        size_t mem_estimate = 232 + (seen_kmers.size() * 8);
        std::cout << "Seen kmer size: " << mem_estimate << std::endl;
        
        local_dict.push_back({seqname, result});
        
        parser.advance();
    }
    
    output_dump(local_dict, out_prefix);
    
    // Note: Saving kmers to binary file not implemented in this version
    // Would require serialization library like Boost.Serialization or custom implementation
}

/**
 * Main deduplication function
 */
void deduplicate(const Args& args) {
    std::cout << "Starting deduplication..." << std::endl;
    
    // Read input files
    std::vector<std::string> fastas;
    
    // Check if input is a list file or direct fasta files
    if (args.input.size() == 1) {
        std::string input = args.input[0];
        if (input.find(".txt") != std::string::npos || input.find(".list") != std::string::npos) {
            // Read list file
            std::ifstream list_file(input);
            std::string line;
            while (std::getline(list_file, line)) {
                if (!line.empty()) {
                    fastas.push_back(line);
                }
            }
        } else {
            fastas = args.input;
        }
    } else {
        fastas = args.input;
    }
    
    // Filter to valid files
    std::vector<std::string> valid_fastas;
    for (const auto& fasta : fastas) {
        if (fs::exists(fasta)) {
            valid_fastas.push_back(fasta);
        } else {
            std::cerr << "Warning: could not find fasta: " << fasta << std::endl;
        }
    }
    
    // Write basename to file map
    std::string basename_fasta_file = args.output_dir + "/basename_fasta_match.txt";
    std::ofstream basename_file(basename_fasta_file);
    for (const auto& fasta : valid_fastas) {
        std::string basename = get_fasta_basename(fasta);
        basename_file << basename << "\t" << fasta << "\n";
    }
    basename_file.close();
    
    // Initialize seen kmers
    std::unordered_set<uint64_t> seen_kmers;
    // Note: Loading from pickle file would require custom deserialization
    
    // Process each fasta
    for (size_t i = 0; i < valid_fastas.size(); ++i) {
        const auto& fasta = valid_fastas[i];
        std::cout << "Deduplicating " << fasta << std::endl;
        
        bool save_kmers = args.save_every > 0 && (i + 1) % args.save_every == 0;
        deduplicate_genome(fasta, seen_kmers, save_kmers, args);
    }
    
    std::cout << "Deduplication complete!" << std::endl;
}

/**
 * Print usage information
 */
void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options] <input_files...>\n\n"
              << "Options:\n"
              << "  -k, --kmer <int>           K-mer size (default: 32)\n"
              << "  -l, --sample-len <int>     Sample length (default: 1000)\n"
              << "  -m, --min-sample-len <int> Minimum sample length (default: same as sample_len)\n"
              << "  -o, --output-dir <path>    Output directory (default: dedup_out/)\n"
              << "  -r, --retain <float>       Likelihood to allow duplicate kmer [0,1] (default: 0.0)\n"
              << "  -s, --save-every <int>     Save seen kmers every n samples (default: 0)\n"
              << "  --no-overlap               Keep neighboring samples discrete\n"
              << "  --no-save-kmers-at-end     Don't save kmers at program end\n"
              << "  --seed <int>               Random seed (default: 123)\n"
              << "  -h, --help                 Show this help message\n";
}

/**
 * Parse command line arguments
 */
Args parse_args(int argc, char* argv[]) {
    Args args;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            exit(0);
        } else if (arg == "-k" || arg == "--kmer") {
            if (i + 1 < argc) {
                args.kmer = std::stoi(argv[++i]);
            }
        } else if (arg == "-l" || arg == "--sample-len") {
            if (i + 1 < argc) {
                args.sample_len = std::stoi(argv[++i]);
            }
        } else if (arg == "-m" || arg == "--min-sample-len") {
            if (i + 1 < argc) {
                args.min_sample_len = std::stoi(argv[++i]);
            }
        } else if (arg == "-o" || arg == "--output-dir") {
            if (i + 1 < argc) {
                args.output_dir = argv[++i];
            }
        } else if (arg == "-r" || arg == "--retain") {
            if (i + 1 < argc) {
                args.retain = std::stod(argv[++i]);
            }
        } else if (arg == "-s" || arg == "--save-every") {
            if (i + 1 < argc) {
                args.save_every = std::stoi(argv[++i]);
            }
        } else if (arg == "--seed") {
            if (i + 1 < argc) {
                args.seed = std::stoi(argv[++i]);
            }
        } else if (arg == "--no-overlap") {
            args.no_overlap = true;
        } else if (arg == "--no-save-kmers-at-end") {
            args.no_save_kmers_at_end = true;
        } else if (arg[0] != '-') {
            args.input.push_back(arg);
        }
    }
    
    // Set min_sample_len to sample_len if not specified
    if (args.min_sample_len == -1) {
        args.min_sample_len = args.sample_len;
    }
    
    return args;
}

/**
 * Main function
 */
int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }
        
        Args args = parse_args(argc, argv);
        
        // Validate arguments
        if (args.input.empty()) {
            std::cerr << "Error: No input files provided" << std::endl;
            return 1;
        }
        
        if (args.kmer < 1) {
            std::cerr << "Error: kmer size must be a positive integer" << std::endl;
            return 1;
        }
        
        if (args.retain < 0.0 || args.retain > 1.0) {
            std::cerr << "Error: retain rate must be between 0.0 and 1.0" << std::endl;
            return 1;
        }
        
        if (args.sample_len < args.kmer) {
            std::cerr << "Warning: sample_len is smaller than k, no deduplication will occur" << std::endl;
        }
        
        if (args.min_sample_len > args.sample_len) {
            std::cerr << "Warning: min_sample_len > sample_len, setting min_sample_len = sample_len" << std::endl;
            args.min_sample_len = args.sample_len;
        }
        
        if (args.kmer < 16) {
            std::cerr << "Warning: small kmer sizes will result in very strict deduplication" << std::endl;
        }
        
        // Check if output directory exists and prompt user
        if (fs::exists(args.output_dir)) {
            std::cout << "Output directory already exists. Continue and potentially overwrite files? (y/n): ";
            std::string response;
            std::getline(std::cin, response);
            if (response != "y" && response != "Y" && response != "yes" && response != "Yes") {
                std::cout << "Exiting..." << std::endl;
                return 1;
            }
        }
        
        // Create output directory if it doesn't exist
        if (!fs::exists(args.output_dir)) {
            fs::create_directories(args.output_dir);
        }
        
        // Initialize random number generator
        rng_gen.seed(args.seed);
        
        // Run deduplication
        deduplicate(args);
        
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
