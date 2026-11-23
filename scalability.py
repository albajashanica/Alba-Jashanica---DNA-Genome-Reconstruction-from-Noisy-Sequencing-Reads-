import sys
import csv
import math
import time
import argparse
import tracemalloc

from simulator import random_genome, simulate_reads, reverse_complement
from indexing import build_oriented_reads, build_kmer_index
from overlap import find_candidates_with_offset, verify_candidates
from assembler import build_graph, greedy_assemble
from consensus import build_pileup, build_consensus, edit_distance_banded
from mutations import mutate_genome, align_to_reference, merge_events, evaluate
from filtering import filter_variants


def reset_peak():
    if hasattr(tracemalloc, "reset_peak"):
        tracemalloc.reset_peak()


def run_size(num_reads, coverage=30, read_len=100, k=12, error_rate=0.01,
             do_consensus=True, var_limit=1200):
    tracemalloc.start()
    genome_len = num_reads * read_len // coverage
    reference = random_genome(genome_len)


    total_mut = genome_len // 300
    if total_mut >= 3:
        n_ins = max(1, total_mut // 4)
        n_del = max(1, total_mut // 4)
        n_sub = max(1, total_mut - n_ins - n_del)
        sample, truth = mutate_genome(reference, n_sub, n_ins, n_del)
    else:
        sample, truth = reference, []

    reads = simulate_reads(sample, num_reads, read_len, error_rate)

    t0 = time.perf_counter()

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    t1 = time.perf_counter()

    candidates = find_candidates_with_offset(oriented, index, k)
    t2 = time.perf_counter()

    overlaps = verify_candidates(oriented, candidates)
    t3 = time.perf_counter()

    graph = build_graph(overlaps)
    contigs = greedy_assemble(oriented, graph)
    t4 = time.perf_counter()

    used, corrected = None, None
    if do_consensus:
        draft = max(contigs, key=len)
        counts, ins, used = build_pileup(draft, reads, k)
        corrected = build_consensus(draft, counts, ins)
    t5 = time.perf_counter()

    peak = tracemalloc.get_traced_memory()[1]
    reset_peak()

    tp = fp = fn = accuracy = ed_sample = None
    n_variants = None
    t_var = None
    mem_var = None
    if do_consensus and num_reads <= var_limit:
        tv0 = time.perf_counter()
        rc = reverse_complement(corrected)
        if edit_distance_banded(rc, reference) < edit_distance_banded(corrected, reference):
            corrected = rc
        events = merge_events(align_to_reference(corrected, reference))
        edge = 15
        detected = [e for e in events if edge <= e[1] <= len(reference) - edge]
        rows = filter_variants(detected, reference, reads)
        final = [r[0] for r in rows if r[3]]
        tp, fp, fn = evaluate(final, truth)
        n_variants = len(final)
        ed_sample = edit_distance_banded(corrected, sample)
        accuracy = round(1 - ed_sample / len(sample), 4)
        t_var = round(time.perf_counter() - tv0, 2)
        mem_var = round(tracemalloc.get_traced_memory()[1] / 1024 / 1024, 1)
    tracemalloc.stop()

    n = len(oriented)
    longest = max(len(c) for c in contigs)
    t_total = t5 - t0
    return {
        "reads": num_reads,
        "genome_len": genome_len,
        "naive": n * (n - 1),
        "candidates": len(candidates),
        "overlaps": len(overlaps),
        "contigs": len(contigs),
        "longest_pct": round(100 * longest / genome_len, 1),
        "total_contig_bases_pct": round(100 * sum(map(len, contigs)) / genome_len, 1),
        "reads_used_pct": round(100 * used / num_reads, 1) if used is not None else None,
        "t_index": round(t1 - t0, 2),
        "t_cand": round(t2 - t1, 2),
        "t_verify": round(t3 - t2, 2),
        "t_graph": round(t4 - t3, 2),
        "t_consensus": round(t5 - t4, 2),
        "t_total": round(t_total, 2),
        "ms_per_read": round(1000 * t_total / num_reads, 1),
        "memory_mb": round(peak / 1024 / 1024, 1),
        "mutations_true": len(truth),
        "variants_found": n_variants,
        "variants_correct": tp,
        "false_alarms": fp,
        "missed": fn,
        "accuracy": accuracy,
        "t_variants": t_var,
        "mem_variants_mb": mem_var,
    }


def fit_power(xs, ys):
    """Pershtat y = a * x^b ne shkalle log-log. Kthen (a, b)."""
    pts = [(math.log(x), math.log(y)) for x, y in zip(xs, ys) if x > 0 and y and y > 0]
    if len(pts) < 2:
        return None
    mx = sum(p[0] for p in pts) / len(pts)
    my = sum(p[1] for p in pts) / len(pts)
    den = sum((p[0] - mx) ** 2 for p in pts)
    if den == 0:
        return None
    b = sum((p[0] - mx) * (p[1] - my) for p in pts) / den
    a = math.exp(my - b * mx)
    return a, b


def analyze_scaling(results):
    reads = [r["reads"] for r in results]
    print("\n=== Analiza e shkallezimit (pershtatje y ~ n^alpha) ===")
    series = [
        ("Krahasime naive", "naive"),
        ("Kandidate me index", "candidates"),
        ("Koha (s)", "t_total"),
        ("Memoria (MB)", "memory_mb"),
    ]
    fits = {}
    for label, key in series:
        fit = fit_power(reads, [r[key] for r in results])
        if fit:
            fits[key] = fit
            print(f"{label:<22} alpha = {fit[1]:.2f}")
    print("(alpha ~ 1 = rritje lineare, alpha ~ 2 = kuadratike)")

    if "t_total" in fits and "memory_mb" in fits:
        print("\nEkstrapolim (jo matje reale, bazuar ne pershtatjen e mesiperm):")
        for target in (500_000, 1_000_000):
            t = fits["t_total"][0] * target ** fits["t_total"][1]
            m = fits["memory_mb"][0] * target ** fits["memory_mb"][1]
            print(f"  {target:>8} reads -> koha ~ {t/60:.0f} min, memoria ~ {m/1024:.1f} GB")


def main():
    ap = argparse.ArgumentParser(description="Vleresim i scalability te rindertimit te genomit")
    ap.add_argument("sizes", nargs="*", type=int,
                    help="numri i reads per secilin test (default: 300 ... 19200)")
    ap.add_argument("--coverage", type=int, default=30)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--k", type=int, default=12)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--var-limit", type=int, default=1200,
                    help="detektimi i variantave ekzekutohet vetem deri ne kete numer reads")
    ap.add_argument("--no-consensus", action="store_true",
                    help="kapercen consensus dhe variant detection")
    ap.add_argument("--out", default="scalability_results.csv")
    args = ap.parse_args()

    sizes = args.sizes or [300, 600, 1200, 2400, 4800, 9600, 19200]
    results = []
    print("reads   genomi   naive        kandidate  overlap   max%    bases%   perdorur%  "
          "koha(s)   ms/read  mem(MB)", flush=True)

    with open(args.out, "w", newline="") as f:
        w = None
        for size in sizes:
            r = run_size(size, args.coverage, args.read_len, args.k, args.error_rate,
                         do_consensus=not args.no_consensus, var_limit=args.var_limit)
            results.append(r)
            if w is None:
                w = csv.DictWriter(f, fieldnames=list(r.keys()))
                w.writeheader()
            w.writerow(r)
            f.flush()

            used = r["reads_used_pct"] if r["reads_used_pct"] is not None else "-"
            print(f"{r['reads']:<7} {r['genome_len']:<8} {r['naive']:<12} {r['candidates']:<10} "
                  f"{r['overlaps']:<9} {r['longest_pct']:<7} {r['total_contig_bases_pct']:<8} {used:<10} "
                  f"{r['t_total']:<9} {r['ms_per_read']:<8} {r['memory_mb']}", flush=True)
            if r["variants_found"] is not None:
                print(f"        variante: {r['variants_correct']}/{r['mutations_true']} te sakta, "
                      f"{r['false_alarms']} alarme te rreme, {r['missed']} te humbura, "
                      f"saktesia={r['accuracy']}, koha variant={r['t_variants']}s", flush=True)

    print("U ruajt:", args.out)

    if len(results) >= 3:
        analyze_scaling(results)


if __name__ == "__main__":
    main()
