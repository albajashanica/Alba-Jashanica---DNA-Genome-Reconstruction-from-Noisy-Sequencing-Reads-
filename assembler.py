from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates

def build_graph(overlaps):
    graph = {}
    for a, b, ov_len, err in overlaps:
        graph.setdefault(a, {})[b] = (ov_len * (1 - err), ov_len)
    return graph

def greedy_assemble(oriented, graph):
    n = len(oriented)
    indeg = [0] * n
    for a, nbrs in graph.items():
        for b in nbrs:
            indeg[b] += 1

    used = [False] * n
    contigs = []
    for start in sorted(range(n), key=lambda x: indeg[x]):
        if used[start]:
            continue
        used[start] = used[start ^ 1] = True
        contig = oriented[start][2]
        cur = start
        while True:
            best, best_w, best_ov = None, -1, 0
            for b, (w, ov_len) in graph.get(cur, {}).items():
                if used[b]:
                    continue
                if len(oriented[b][2]) - ov_len < 1:   # b brenda cur, s'zgjeron
                    continue
                if w > best_w:
                    best, best_w, best_ov = b, w, ov_len
            if best is None:
                break
            contig += oriented[best][2][best_ov:]
            used[best] = used[best ^ 1] = True
            cur = best
        contigs.append(contig)
    return contigs


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