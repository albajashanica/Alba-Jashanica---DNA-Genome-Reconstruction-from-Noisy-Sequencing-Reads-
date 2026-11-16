import argparse
import time
import tracemalloc

from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble
from consensus import build_pileup, build_consensus, edit_distance_banded
from mutations import mutate_genome, align_to_reference, merge_events, evaluate
from filtering import filter_variants

def run_pipeline(genome_len=2000, num_reads=300, read_len=100, error_rate=0.01,
                 k=12, min_overlap=30, max_overlap_error=0.08):
    
    reference = random_genome(genome_len)
    sample, truth = mutate_genome(reference)
    reads = simulate_reads(sample, num_reads, read_len, error_rate)
    print("reads:", len(reads)) 

if __name__ == "__main__":  
    run_pipeline()  