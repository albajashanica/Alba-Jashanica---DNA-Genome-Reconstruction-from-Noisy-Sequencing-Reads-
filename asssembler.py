from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates

def build_graph(overlaps):
    graph = {}
    for a, b, ov_len, err in overlaps:
        graph.setdefault(a, {})[b] = (ov_len * (1 - err), ov_len)
    return graph

