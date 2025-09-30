import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dedup_functions import *

test_string = "ACGTACGTACGT"
test_k = 3
test_sample_len = 5
test_min_sample_len = 4
test_seen_kmers = set()


def test_encode_kmer():
    assert encode_kmer("ACG") == 6
    assert encode_kmer("TCG") == 54
    assert encode_kmer("AAA") == 0
    assert encode_kmer("AAC") != encode_kmer("ACA")

def test_decode_kmer():
    assert decode_kmer(6, 3) == "ACG"
    assert decode_kmer(54, 3) == "TCG"
    assert decode_kmer(0, 3) == "AAA"
    assert decode_kmer(encode_kmer("ACA"), 3) == "ACA"

def test_n_check():
    assert n_check("ACGTCGGT", 0, 7, []) == ([], 0, 7)
    assert n_check("ACGTNGGT", 0, 7, []) == ([[0,5]], 5, 7)
    assert n_check("ANGTNGGT", 0, 7, []) == ([[0,5]], 5, 7)


test_encode_kmer()
test_decode_kmer()
test_n_check()