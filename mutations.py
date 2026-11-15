import random
from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble
from consensus import build_pileup, build_consensus, edit_distance_banded

def mutate_genome(genome, n_sub=5, n_ins=3, n_del=3, seed=7):
    rng = random.Random(seed)
    total = n_sub + n_ins + n_del
    positions = []
    while len(positions) < total:
        p = rng.randint(60, len(genome) - 60)
        if all(abs(p - q) >= 60 for q in positions):
            positions.append(p)
    kinds = ["SUB"] * n_sub + ["INS"] * n_ins + ["DEL"] * n_del
    rng.shuffle(kinds)
    truth = []
    for p, kind in zip(positions, kinds):
        if kind == "SUB":
            alt = rng.choice([b for b in "ACGT" if b != genome[p]])
            truth.append(("SUB", p, genome[p], alt))
        elif kind == "INS":
            truth.append(("INS", p, "", rng.choice("ACGT")))
        else:
            truth.append(("DEL", p, genome[p], ""))
    sample = genome
    for kind, p, ref, alt in sorted(truth, key=lambda t: -t[1]):
        if kind == "SUB":
            sample = sample[:p] + alt + sample[p + 1:]
        elif kind == "INS":
            sample = sample[:p] + alt + sample[p:]
        else:
            sample = sample[:p] + sample[p + 1:]
    return sample, sorted(truth, key=lambda t: t[1])

def align_to_reference(cons, ref, band=60):
    n, m = len(cons), len(ref)
    INF = 10 ** 9
    dp = {(0, 0): 0}
    for j in range(1, min(m, band) + 1):
        dp[(0, j)] = j
    for i in range(1, n + 1):
        for j in range(max(0, i - band), min(m, i + band) + 1):
            best = dp.get((i - 1, j), INF) + 1
            if j > 0:
                cost = 0 if cons[i - 1] == ref[j - 1] else 1
                best = min(best,
                           dp.get((i - 1, j - 1), INF) + cost,
                           dp.get((i, j - 1), INF) + 1)
            dp[(i, j)] = best
    i, j = n, m
    events = []
    while i > 0 or j > 0:
        cur = dp[(i, j)]
        if i > 0 and j > 0:
            cost = 0 if cons[i - 1] == ref[j - 1] else 1
            if cur == dp.get((i - 1, j - 1), INF) + cost:
                if cost:
                    events.append(("SUB", j - 1, ref[j - 1], cons[i - 1]))
                i -= 1
                j -= 1
                continue
        if i > 0 and cur == dp.get((i - 1, j), INF) + 1:
            events.append(("INS", j, "", cons[i - 1]))
            i -= 1
        else:
            events.append(("DEL", j - 1, ref[j - 1], ""))
            j -= 1
    events.reverse()
    return events
def merge_events(events):
    merged = []
    for kind, p, ref, alt in events:
        if merged:
            lk, lp, lref, lalt = merged[-1]
            if kind == "DEL" and lk == "DEL" and p == lp + len(lref):
                merged[-1] = (lk, lp, lref + ref, lalt)
                continue
            if kind == "INS" and lk == "INS" and p == lp:
                merged[-1] = (lk, lp, lref, lalt + alt)
                continue
        merged.append((kind, p, ref, alt))
    return merged


def evaluate(detected, truth, tol=3):
    matched = set()
    tp = 0
    for d in detected:
        for idx, t in enumerate(truth):
            if idx in matched:
                continue
            if d[0] == t[0] and abs(d[1] - t[1]) <= tol and (d[0] != "SUB" or d[3] == t[3]):
                matched.add(idx)
                tp += 1
                break
    return tp, len(detected) - tp, len(truth) - tp

if __name__ == "__main__":
    k = 12
    reference = random_genome(2000)
    sample, truth = mutate_genome(reference)
    reads = simulate_reads(sample, num_reads=300)

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates)
    graph = build_graph(overlaps)
    contigs = sorted(greedy_assemble(oriented, graph), key=len, reverse=True)

    counts, ins, used = build_pileup(contigs[0], reads, k)
    corrected = build_consensus(contigs[0], counts, ins)

    rc = reverse_complement(corrected)
    if edit_distance_banded(rc, reference) < edit_distance_banded(corrected, reference):
        corrected = rc

    events = merge_events(align_to_reference(corrected, reference))
    edge = 15
    detected = [e for e in events if edge <= e[1] <= len(reference) - edge]
    tp, fp, fn = evaluate(detected, truth)