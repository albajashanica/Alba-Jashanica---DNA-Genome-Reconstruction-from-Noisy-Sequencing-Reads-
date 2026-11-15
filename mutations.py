
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

if __name__ == "__main__":
    k = 12
    reference = random_genome(2000)
    sample, truth = mutate_genome(reference)
    reads = simulate_reads(sample, num_reads=300)
