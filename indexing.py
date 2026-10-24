from collections import defaultdict
from simulator import random_genome, simulate_reads, reverse_complement

def build_oriented_reads(reads):
    oriented = []
    for i, r in enumerate(reads):
        oriented.append((i, "+", r))
        oriented.append((i, "-", reverse_complement(r)))
    return oriented

def build_kmer_index(oriented, k):
    index = defaultdict(list)
    for node_id, (_, _, seq) in enumerate(oriented):
        for pos in range(len(seq) - k + 1):
            index[seq[pos:pos + k]].append((node_id, pos))
    return index