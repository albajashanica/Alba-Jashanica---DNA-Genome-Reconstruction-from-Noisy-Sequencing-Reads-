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

def banded_overlap(a_suffix, b, band=5):
    m, n = len(a_suffix), len(b)
    INF = 10**9
    prev = [INF] * (n + 1)
    for j in range(0, min(band, n) + 1):
        prev[j] = j
    for i in range(1, m + 1):
        cur = [INF] * (n + 1)
        lo, hi = max(0, i - band), min(n, i + band)
        for j in range(lo, hi + 1):
            best = prev[j] + 1
            if j > 0:
                cost = 0 if a_suffix[i - 1] == b[j - 1] else 1
                best = min(best, prev[j - 1] + cost, cur[j - 1] + 1)
            cur[j] = best
        prev = cur
    lo, hi = max(0, m - band), min(n, m + band)
    if lo > hi:
        return None
    j_best = min(range(lo, hi + 1), key=lambda j: prev[j])
    return prev[j_best]


def verify_candidates(oriented, candidates, min_overlap=30, max_error=0.08):
    overlaps = []
    for (a, b), off in candidates.items():
        seq_a, seq_b = oriented[a][2], oriented[b][2]
        ov_len = len(seq_a) - off
        if ov_len < min_overlap:
            continue
        edits = banded_overlap(seq_a[off:], seq_b)
        if edits is None:
            continue
        error_rate = edits / ov_len
        if error_rate <= max_error:
            overlaps.append((a, b, ov_len, error_rate))
    return overlaps

if __name__ == "__main__":
    k = 12
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates)

    print("Kandidate:", len(candidates))
    print("Overlap te pranuar:", len(overlaps))
    if overlaps:
        avg_len = sum(o[2] for o in overlaps) / len(overlaps)
        avg_q = sum(1 - o[3] for o in overlaps) / len(overlaps)
        print("Gjatesia mesatare e overlap:", round(avg_len, 1))
        print("Cilesia mesatare (1 - error rate):", round(avg_q, 3))