from collections import defaultdict
from simulator import random_genome, simulate_reads
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates

def build_overlap_graph(overlaps):
    graph = defaultdict(dict)   # graph[a][b] = pesha
    for a, b, ov_len, err in overlaps:
        graph[a][b] = ov_len * (1 - err)
    return graph

if __name__ == "__main__":
    k = 12
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates)

    graph = build_overlap_graph(overlaps)
    n = len(oriented)
    
    print("Nyje:", n)
    print("Grafi u ndërtua me sukses!")