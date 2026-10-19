import random

def random_genome(length, seed=42):
    random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(length))

def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def add_errors(read, error_rate):
    result = []
    for base in read:
        r = random.random()
        if r < error_rate / 3:            
            result.append(random.choice([b for b in "ACGT" if b != base]))
        elif r < 2 * error_rate / 3:      
            result.append(base)
            result.append(random.choice("ACGT"))
        elif r < error_rate:              
            continue
        else:
            result.append(base)
    return "".join(result)