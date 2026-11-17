import sys
import csv
import time
import tracemalloc

from simulator import random_genome, simulate_reads

def run_size(num_reads, coverage=15, read_len=100, k=12):
    genome_len = num_reads * read_len // coverage
    genome = random_genome(genome_len)
    reads = simulate_reads(genome, num_reads, read_len, 0.01)
    
if __name__ == "__main__":
    print("Simulimi u realizua me sukses.")