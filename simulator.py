import random

def random_genome(length, seed=42):
    random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(length))

def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

