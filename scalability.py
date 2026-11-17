import sys
import csv
import time
import tracemalloc

from simulator import random_genome, simulate_reads
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates

def run_size(num_reads, coverage=15, read_len=100, k=12):
    genome_len = num_reads * read_len // coverage
    genome = random_genome(genome_len)
    reads = simulate_reads(genome, num_reads, read_len, 0.01)
    
    tracemalloc.start()
    t0 = time.perf_counter()

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    t1 = time.perf_counter()

    candidates = find_candidates_with_offset(oriented, index, k)
    t2 = time.perf_counter()

    overlaps = verify_candidates(oriented, candidates)
    t3 = time.perf_counter()

    tracemalloc.stop()

    return {
        "reads": num_reads,
        "genome_len": genome_len,
        "naive": len(oriented) * (len(oriented) - 1),
        "candidates": len(candidates),
        "overlaps": len(overlaps),
    }
if __name__ == "__main__":
    print("Simulimi u realizua me sukses.")
    r = run_size(300)
    print(r)