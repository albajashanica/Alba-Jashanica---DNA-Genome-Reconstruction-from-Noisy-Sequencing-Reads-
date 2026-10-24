from collections import defaultdict
from simulator import random_genome, simulate_reads
from indexing import build_oriented_reads, build_kmer_index

def find_candidates_with_offset(oriented, index, k, min_shared=3):
    candidates = {}
    for a, (read_a, _, seq_a) in enumerate(oriented):
        votes = defaultdict(int)
        for pos_a in range(len(seq_a) - k + 1):
            for b, pos_b in index[seq_a[pos_a:pos_a + k]]:
                if oriented[b][0] != read_a and pos_a >= pos_b:
                    votes[(b, pos_a - pos_b)] += 1
        best = {}
        for (b, off), c in votes.items():
            if b not in best or c > best[b][1]:
                best[b] = (off, c)
        for b, (off, c) in best.items():
            if c >= min_shared:
                candidates[(a, b)] = off
    return candidates

if __name__ == "__main__":
    k = 12
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)
    
    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    
    print("Kandidate:", len(candidates))