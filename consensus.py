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
            ops.append(("M", j - 1, read[i - 1], i - 1)); i -= 1; j -= 1
        elif dp[i][j] == dp[i - 1][j] + 1:
            ops.append(("I", j - 1, read[i - 1], i - 1)); i -= 1
        else:
            ops.append(("D", j - 1, "-", None)); j -= 1
    return edits, ops

def align_read(read, contig, cindex, k, band=10, quality=None):
    best = None
    for seq, qual in ((read, quality), (reverse_complement(read), quality[::-1] if quality else None)):
        votes = defaultdict(int)
        for pr in range(len(seq) - k + 1):
            for pc in cindex.get(seq[pr:pr + k], []):
                votes[pc - pr] += 1
        if votes:
            off, v = max(votes.items(), key=lambda x: x[1])
            if best is None or v > best[0]:
                best = (v, off, seq, qual)
    if best is None or best[0] < 3:
        return None
    _, off, seq, qual = best
    start = max(0, off - band)
    end = min(len(contig), off + len(seq) + band)
    if end - start < len(seq) // 2:
        return None
    edits, ops = semi_global(seq, contig[start:end])
    if edits > 0.15 * len(seq):
        return None
    return start, ops, qual

def build_pileup(contig, reads, k, qualities=None):
    cindex = defaultdict(list)
    for p in range(len(contig) - k + 1):
        cindex[contig[p:p + k]].append(p)
    counts = [defaultdict(int) for _ in contig]
    ins = defaultdict(lambda: defaultdict(int))
    used = 0
    for ri, read in enumerate(reads):
        qual = qualities[ri] if qualities is not None else None
        res = align_read(read, contig, cindex, k, quality=qual)
        if res is None:
            continue
        used += 1
        start, ops, qual = res
        for t, jj, base, read_pos in ops:
            pos = start + jj
            weight = 1.0
            if qual is not None and read_pos is not None:
                q = max(0, ord(qual[read_pos]) - 33)
                weight = max(0.01, 1 - 10 ** (-q / 10))
            if t == "I":
                ins[pos][base] += weight
            elif 0 <= pos < len(contig):
                counts[pos][base] += weight
    return counts, ins, used

def build_consensus(contig, counts, ins):
    out = []
    for p, base in enumerate(contig):
        c = counts[p]
        total = sum(c.values())
        if total == 0:
            out.append(base)
            continue
        best = max(c, key=c.get)
        if best != "-":
            out.append(best)
        if p in ins:
            b, cnt = max(ins[p].items(), key=lambda x: x[1])
            if cnt > total / 2:
                out.append(b)
    return "".join(out)

def edit_distance_banded(a, b, band=60):
    n, m = len(a), len(b)
    band = max(band, abs(n - m))
    INF = 10 ** 9
    prev = {j: j for j in range(0, min(m, band) + 1)}
    for i in range(1, n + 1):
        cur = {}
        for j in range(max(0, i - band), min(m, i + band) + 1):
            best = prev.get(j, INF) + 1
            if j > 0:
                best = min(best,
                           prev.get(j - 1, INF) + (a[i - 1] != b[j - 1]),
                           cur.get(j - 1, INF) + 1)
            cur[j] = best
        prev = cur
    return prev.get(m, INF)

def distance_to_genome(seq, genome):
    return min(edit_distance_banded(seq, genome),
               edit_distance_banded(reverse_complement(seq), genome))

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

    counts, ins, used = build_pileup(draft, reads, k)
    corrected = build_consensus(draft, counts, ins)

    print("Reads te shtrira mbi contig:", used, "nga", len(reads))
    print("Mbulimi mesatar:", round(sum(sum(c.values()) for c in counts) / len(draft), 1))
    print("Gjatesia draft / e korrigjuar:", len(draft), "/", len(corrected))
    d1 = distance_to_genome(draft, genome)
    d2 = distance_to_genome(corrected, genome)
    print("Edit distance para korrigjimit:", d1)
    print("Edit distance pas korrigjimit:", d2)
    print("Saktesia pas korrigjimit:", round(1 - d2 / len(genome), 4))
