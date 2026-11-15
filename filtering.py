from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble
from consensus import build_pileup, build_consensus, edit_distance_banded
from mutations import mutate_genome, align_to_reference, merge_events, evaluate

def variant_support(variant, reference, reads, flank=12):
    kind, p, ref_allele, alt_allele = variant
    end = p + len(ref_allele)
    ref_win = reference[p - flank:end + flank]
    alt_win = reference[p - flank:p] + alt_allele + reference[end:end + flank]
    ref_support = alt_support = 0
    for r in reads:
        both = (r, reverse_complement(r))
        if any(ref_win in s for s in both):
            ref_support += 1
        if any(alt_win in s for s in both):
            alt_support += 1
    return alt_support, ref_support


def filter_variants(variants, reference, reads, min_alt=3):
    rows = []
    for v in variants:
        alt_s, ref_s = variant_support(v, reference, reads)
        passed = alt_s >= min_alt and alt_s > ref_s
        rows.append((v, alt_s, ref_s, passed))
    return rows


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

    rows = filter_variants(detected, reference, reads)
    print("Variant | reads me variantin | reads me referencen | rezultati")
    for v, a, r, ok in rows:
        print(v, "|", a, "|", r, "|", "PASS" if ok else "FILTERED")

    final = [v for v, a, r, ok in rows if ok]
