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
