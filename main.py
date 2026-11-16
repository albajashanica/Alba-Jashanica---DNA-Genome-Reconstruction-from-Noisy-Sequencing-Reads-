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
    
    
    tracemalloc.start()
    t0 = time.perf_counter()

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates, min_overlap, max_overlap_error)
    graph = build_graph(overlaps)
    contigs = sorted(greedy_assemble(oriented, graph), key=len, reverse=True)
    draft = contigs[0]

    counts, ins, used = build_pileup(draft, reads, k)
    corrected = build_consensus(draft, counts, ins)

    rc = reverse_complement(corrected)
    if edit_distance_banded(rc, reference) < edit_distance_banded(corrected, reference):
        corrected = rc
    events = merge_events(align_to_reference(corrected, reference))
    edge = 15
    detected = [e for e in events if edge <= e[1] <= len(reference) - edge]
    rows = filter_variants(detected, reference, reads)
    final_rows = [r for r in rows if r[3]]
    final = [r[0] for r in final_rows]

    elapsed = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    return corrected, final_rows

if __name__ == "__main__":  
    run_pipeline()  