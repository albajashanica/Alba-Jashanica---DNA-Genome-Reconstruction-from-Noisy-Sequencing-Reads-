from collections import defaultdict
from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble

def semi_global(read, window):
    m, w = len(read), len(window)
    dp = [[0] * (w + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i
        for j in range(1, w + 1):
            cost = 0 if read[i - 1] == window[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j - 1] + cost, dp[i - 1][j] + 1, dp[i][j - 1] + 1)
    j = min(range(w + 1), key=lambda x: dp[m][x])
    edits = dp[m][j]
    i = m
    ops = []  
    while i > 0:
        if j > 0 and dp[i][j] == dp[i - 1][j - 1] + (0 if read[i - 1] == window[j - 1] else 1):
            ops.append(("M", j - 1, read[i - 1])); i -= 1; j -= 1
        elif dp[i][j] == dp[i - 1][j] + 1:
            ops.append(("I", j - 1, read[i - 1])); i -= 1
        else:
            ops.append(("D", j - 1, "-")); j -= 1
    return edits, ops

def align_read(read, contig, cindex, k, band=10):
    best = None
    for seq in (read, reverse_complement(read)):
        votes = defaultdict(int)
        for pr in range(len(seq) - k + 1):
            for pc in cindex.get(seq[pr:pr + k], []):
                votes[pc - pr] += 1
        if votes:
            off, v = max(votes.items(), key=lambda x: x[1])
            if best is None or v > best[0]:
                best = (v, off, seq)
    if best is None or best[0] < 3:
        return None
    _, off, seq = best
    start = max(0, off - band)
    end = min(len(contig), off + len(seq) + band)
    if end - start < len(seq) // 2:
        return None
    edits, ops = semi_global(seq, contig[start:end])
    if edits > 0.15 * len(seq):
        return None
    return start, ops

if __name__ == "__main__":
    k = 12
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates)
    graph = build_graph(overlaps)
    contigs = sorted(greedy_assemble(oriented, graph), key=len, reverse=True)
    draft = contigs[0]