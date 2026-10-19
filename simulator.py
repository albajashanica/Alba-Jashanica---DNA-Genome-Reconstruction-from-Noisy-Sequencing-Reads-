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

def simulate_reads(genome, num_reads, read_len=100, error_rate=0.01, rc_prob=0.5):
    reads = []
    for _ in range(num_reads):
        start = random.randint(0, len(genome) - read_len)
        read = genome[start:start + read_len]
        read = add_errors(read, error_rate)
        if random.random() < rc_prob:
            read = reverse_complement(read)
        reads.append(read)
    return reads

if __name__ == "__main__":
    genome = random_genome(2000)
    reads = simulate_reads(genome, num_reads=300)
    print("Gjatesia e genomit:", len(genome))
    print("Numri i reads:", len(reads))
    print("Read i pare:", reads[0])