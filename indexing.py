from collections import defaultdict
from simulator import random_genome, simulate_reads, reverse_complement

def build_oriented_reads(reads):
    oriented = []
    for i, r in enumerate(reads):
        oriented.append((i, "+", r))
        oriented.append((i, "-", reverse_complement(r)))
    return oriented

def build_kmer_index(oriented, k, max_occurrences=200):
    """Build an inverted index, omitting highly repetitive k-mers.

    Common k-mers create enormous candidate lists and usually occur in repeats
    that cannot uniquely place a read. Skipping those buckets bounds candidate
    generation; callers can raise max_occurrences when repeat sensitivity is
    more important than speed.
    """
    if k < 1:
        raise ValueError("k must be positive")
    index = defaultdict(list)
    for node_id, (_, _, seq) in enumerate(oriented):
        for pos in range(len(seq) - k + 1):
            index[seq[pos:pos + k]].append((node_id, pos))
    return {mer: hits for mer, hits in index.items() if len(hits) <= max_occurrences}

def find_candidates(oriented, index, k, min_shared=3):
    candidates = {}
    for a, (read_a, _, seq_a) in enumerate(oriented):
        shared = defaultdict(int)
        for pos_a in range(len(seq_a) - k + 1):
            for b, pos_b in index[seq_a[pos_a:pos_a + k]]:
                if oriented[b][0] != read_a and pos_a >= pos_b:
                    shared[b] += 1
        for b, count in shared.items():
            if count >= min_shared:
                candidates[(a, b)] = count
    return candidates

if __name__ == "__main__":
    k = 12
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates(oriented, index, k)

    n = len(oriented)
    print("Nyje (reads x 2 orientime):", n)
    print("K-mers unike ne index:", len(index))
    print("Krahasime ne menyren naive:", n * (n - 1))
    print("Kandidate me index:", len(candidates))
