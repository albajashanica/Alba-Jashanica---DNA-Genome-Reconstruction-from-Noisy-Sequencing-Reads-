import random

def random_genome(length, seed=42):
    random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(length))