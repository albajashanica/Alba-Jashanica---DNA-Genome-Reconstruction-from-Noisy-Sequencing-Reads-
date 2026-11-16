import argparse
import time
import tracemalloc

from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble
from consensus import build_pileup, build_consensus, edit_distance_banded
from mutations import mutate_genome, align_to_reference, merge_events, evaluate
from filtering import filter_variants

def run_pipeline(genome_len=2000, num_reads=300, read_len=100, error_rate=0.01,
                 k=12, min_overlap=30, max_overlap_error=0.08):
    
    reference = random_genome(genome_len)
    sample, truth = mutate_genome(reference)
    reads = simulate_reads(sample, num_reads, read_len, error_rate)


    tracemalloc.start()
    t0 = time.perf_counter()

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates, min_overlap, max_overlap_error)
    graph = build_graph(overlaps)
    contigs = sorted(greedy_assemble(oriented, graph), key=len, reverse=True)
    draft = contigs[0]

    counts, ins, used = build_pileup(draft, reads, k)
    corrected = build_consensus(draft, counts, ins)

    rc = reverse_complement(corrected)
    if edit_distance_banded(rc, reference) < edit_distance_banded(corrected, reference):
        corrected = rc
    events = merge_events(align_to_reference(corrected, reference))
    edge = 15
    detected = [e for e in events if edge <= e[1] <= len(reference) - edge]
    rows = filter_variants(detected, reference, reads)
    final_rows = [r for r in rows if r[3]]
    final = [r[0] for r in final_rows]

    elapsed = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    
    n = len(oriented)
    ed_sample = edit_distance_banded(corrected, sample)
    ed_ref = edit_distance_banded(corrected, reference)
    tp, fp, fn = evaluate(final, truth)
    stats = {
        "reads_total": num_reads,
        "reads_used": used,
        "pct_reads_used": round(100 * used / num_reads, 1),
        "nodes": n,
        "naive_comparisons": n * (n - 1),
        "candidate_comparisons": len(candidates),
        "overlaps": len(overlaps),
        "avg_overlap_len": round(sum(o[2] for o in overlaps) / max(1, len(overlaps)), 1),
        "avg_overlap_quality": round(sum(1 - o[3] for o in overlaps) / max(1, len(overlaps)), 3),
        "contigs": len(contigs),
        "reconstruction_length": len(corrected),
        "true_length": len(sample),
        "edit_distance_sample": ed_sample,
        "accuracy": round(1 - ed_sample / len(sample), 4),
        "edit_distance_reference": ed_ref,
        "mutations_inserted": len(truth),
        "mutations_correct": tp,
        "false_alarms": fp,
        "missed": fn,
        "time_s": round(elapsed, 2),
        "peak_memory_mb": round(peak / 1024 / 1024, 2),
    }
    return corrected, final_rows, stats

def save_outputs(corrected, final_rows):
    with open("reconstructed.fasta", "w") as f:
        f.write(">reconstructed_genome length=" + str(len(corrected)) + "\n")
        for i in range(0, len(corrected), 70):
            f.write(corrected[i:i + 70] + "\n")
    with open("variants_report.txt", "w") as f:
        f.write("lloji\tpozicioni\tbaza_ref\tbaza_re\treads_variant\treads_ref\n")
        for (kind, pos, ref, alt), a, r, _ in final_rows:
            f.write(f"{kind}\t{pos}\t{ref or '-'}\t{alt or '-'}\t{a}\t{r}\n")

def print_report(s):
    print("====== RAPORTI I RINDERTIMIT ======")
    print(f"Reads te perdorura:            {s['reads_used']} / {s['reads_total']} ({s['pct_reads_used']}%)")
    print(f"Nyje ne graf (2 orientime):    {s['nodes']}")
    print(f"Krahasime naive (n^2):         {s['naive_comparisons']}")
    print(f"Krahasime kandidate (index):   {s['candidate_comparisons']}")
    print(f"Overlap te identifikuara:      {s['overlaps']}")
    print(f"Gjatesia mesatare e overlap:   {s['avg_overlap_len']}")
    print(f"Cilesia mesatare e overlap:    {s['avg_overlap_quality']}")
    print(f"Numri i contigs (draft):       {s['contigs']}")
    print(f"Gjatesia e rindertimit:        {s['reconstruction_length']} (e vertete: {s['true_length']})")
    print(f"Edit distance vs mostra:       {s['edit_distance_sample']}")
    print(f"Saktesia vs mostra:            {s['accuracy']}")
    print(f"Edit distance vs referenca:    {s['edit_distance_reference']}")
    print(f"Mutacione te futura:           {s['mutations_inserted']}")
    print(f"Mutacione te zbuluara sakte:   {s['mutations_correct']}")
    print(f"Alarme te rreme / te humbura:  {s['false_alarms']} / {s['missed']}")
    print(f"Koha e ekzekutimit:            {s['time_s']} s")
    print(f"Memoria maksimale (Python):    {s['peak_memory_mb']} MB")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="DNA genome reconstruction")
    ap.add_argument("--genome-len", type=int, default=2000)
    ap.add_argument("--reads", type=int, default=300)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--k", type=int, default=12)
    ap.add_argument("--min-overlap", type=int, default=30)
    ap.add_argument("--max-error", type=float, default=0.08)
    args = ap.parse_args()

    corrected, final_rows, stats = run_pipeline(
        args.genome_len, args.reads, args.read_len, args.error_rate,
        args.k, args.min_overlap, args.max_error)
    print_report(stats)
    save_outputs(corrected, final_rows)
    print("U ruajten: reconstructed.fasta dhe variants_report.txt")