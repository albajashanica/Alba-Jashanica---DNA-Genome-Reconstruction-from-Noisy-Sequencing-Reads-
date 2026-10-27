from collections import defaultdict
from simulator import random_genome, simulate_reads
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates

def build_overlap_graph(overlaps):
    graph = defaultdict(dict)   # graph[a][b] = pesha
    for a, b, ov_len, err in overlaps:
        graph[a][b] = ov_len * (1 - err)
    return graph

def graph_stats(graph, num_nodes):
    num_edges = sum(len(nbrs) for nbrs in graph.values())
    no_out = sum(1 for n in range(num_nodes) if not graph.get(n))
    mirrored = sum(
        1 for a, nbrs in graph.items() for b in nbrs
        if (a ^ 1) in graph.get(b ^ 1, {})
    )
    return num_edges, no_out, mirrored

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
    edges, no_out, mirrored = graph_stats(graph, n)

    print("Skaje:", edges)
    print("Shkalla mesatare dalese:", round(edges / n, 2))
    print("Nyje pa skaj dales:", no_out)
    print("Skaje me pasqyre:", mirrored, "nga", edges)