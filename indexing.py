from collections import defaultdict
from simulator import random_genome, simulate_reads, reverse_complement

def build_oriented_reads(reads):
    oriented = []
    for i, r in enumerate(reads):
        oriented.append((i, "+", r))
        oriented.append((i, "-", reverse_complement(r)))
    return oriented