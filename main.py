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


def read_sequences(path, return_qualities=False):
    """Read DNA sequences from FASTA or FASTQ (plain or .gz)."""
    import gzip
    opener = gzip.open if path.lower().endswith(".gz") else open
    seqs, qualities = [], []
    with opener(path, "rt", encoding="ascii") as fh:
        first = next((line.strip() for line in fh if line.strip()), "")
        if not first:
            raise ValueError(f"No sequences found in {path}")
        if first.startswith(">"):
            chunks = []
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if chunks:
                        seqs.append("".join(chunks).upper())
                        qualities.append(None)
                        chunks = []
                else:
                    chunks.append(line)
            if chunks:
                seqs.append("".join(chunks).upper())
                qualities.append(None)
        elif first.startswith("@"):
            header = first
            record = 1
            while header:
                sequence = fh.readline().strip()
                plus = fh.readline().strip()
                quality = fh.readline().rstrip("\r\n")
                if not sequence or not plus.startswith("+") or len(sequence) != len(quality):
                    raise ValueError(f"Invalid FASTQ record {record}")
                seqs.append(sequence.upper())
                qualities.append(quality)
                header = next((line.strip() for line in fh if line.strip()), "")
                if header and not header.startswith("@"):
                    raise ValueError(f"Invalid FASTQ header before record {record + 1}")
                record += 1
        else:
            raise ValueError("Input must be FASTA or FASTQ (optionally .gz)")
    valid = set("ACGT")
    kept = [(s, q) for s, q in zip(seqs, qualities) if s and set(s) <= valid]
    seqs = [s for s, _ in kept]
    qualities = [q for _, q in kept]
    if not seqs:
        raise ValueError("No valid DNA reads (A/C/G/T) found")
    return (seqs, qualities) if return_qualities else seqs


def run_pipeline(genome_len=333333, num_reads=100000, read_len=100, error_rate=0.01,
                 k=12, min_overlap=30, max_overlap_error=0.08,
                 reads_path=None, reference_path=None):
    if reference_path and not reads_path:
        raise ValueError("--reference is used with --reads-file")
    if not 0 <= max_overlap_error <= 1:
        raise ValueError("max_overlap_error must be between 0 and 1")
    tracemalloc.start()
    t0 = time.perf_counter()
    simulated = reads_path is None
    truth = []
    sample = None
    if simulated:
        reference = random_genome(genome_len)
        sample, truth = mutate_genome(reference)
        reads = simulate_reads(sample, num_reads, read_len, error_rate)
        qualities = None
    else:
        reads, qualities = read_sequences(reads_path, return_qualities=True)
        reference = read_sequences(reference_path)[0] if reference_path else None

    oriented = build_oriented_reads(reads)
    index = build_kmer_index(oriented, k)
    candidates = find_candidates_with_offset(oriented, index, k)
    overlaps = verify_candidates(oriented, candidates, min_overlap, max_overlap_error, k=k)
    graph = build_graph(overlaps)
    contigs = sorted(greedy_assemble(oriented, graph), key=len, reverse=True)
    if not contigs:
        raise ValueError("No contigs were assembled; check read lengths and parameters")
    draft = contigs[0]
    counts, ins, used = build_pileup(draft, reads, k, qualities)
    corrected = build_consensus(draft, counts, ins)

    final_rows = []
    ed_sample = ed_ref = accuracy = accuracy_ref = None
    tp = fp = fn = None
    if reference is not None:
        rc = reverse_complement(corrected)
        if edit_distance_banded(rc, reference) < edit_distance_banded(corrected, reference):
            corrected = rc
        events = merge_events(align_to_reference(corrected, reference))
        edge = min(15, max(0, len(reference) // 10))
        detected = [e for e in events if edge <= e[1] < len(reference) - edge]
        rows = filter_variants(detected, reference, reads)
        final_rows = [r for r in rows if r[3]]
        ed_ref = edit_distance_banded(corrected, reference)
        accuracy_ref = round(1 - ed_ref / max(1, len(reference)), 4)
        if simulated:
            ed_sample = edit_distance_banded(corrected, sample)
            accuracy = round(1 - ed_sample / max(1, len(sample)), 4)
            tp, fp, fn = evaluate([r[0] for r in final_rows], truth)

    elapsed = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    n = len(oriented)
    stats = {
        "reads_total": len(reads), "reads_used": used,
        "pct_reads_used": round(100 * used / max(1, len(reads)), 1),
        "nodes": n, "naive_comparisons": n * (n - 1),
        "candidate_comparisons": len(candidates), "overlaps": len(overlaps),
        "avg_overlap_len": round(sum(o[2] for o in overlaps) / max(1, len(overlaps)), 1),
        "avg_overlap_quality": round(sum(1 - o[3] for o in overlaps) / max(1, len(overlaps)), 3),
        "contigs": len(contigs), "reconstruction_length": len(corrected),
        "assembled_bases": sum(map(len, contigs)) - len(draft) + len(corrected),
        "true_length": len(sample) if sample is not None else (len(reference) if reference else None),
        "edit_distance_sample": ed_sample, "accuracy": accuracy,
        "edit_distance_reference": ed_ref, "accuracy_reference": accuracy_ref,
        "mutations_inserted": len(truth) if simulated else None,
        "mutations_correct": tp, "false_alarms": fp, "missed": fn,
        "time_s": round(elapsed, 2), "peak_memory_mb": round(peak / 1024 / 1024, 2),
        "reference_available": reference is not None,
    }
    output_contigs = [corrected] + contigs[1:]
    return corrected, final_rows, stats, output_contigs


def save_outputs(corrected, final_rows, contigs=None):
    contigs = contigs or [corrected]
    with open("reconstructed.fasta", "w", encoding="ascii") as f:
        for ci, sequence in enumerate(contigs, 1):
            f.write(f">contig_{ci} length={len(sequence)}\n")
            for i in range(0, len(sequence), 70):
                f.write(sequence[i:i + 70] + "\n")
    with open("variants_report.txt", "w", encoding="utf-8") as f:
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
    print(f"Gjatesia e rindertimit:        {s['reconstruction_length']} (referenca: {s['true_length']})")
    print(f"Baza gjithsej ne contig-e:     {s['assembled_bases']}")
    print(f"Edit distance vs mostra:       {s['edit_distance_sample']}")
    print(f"Saktesia vs mostra:            {s['accuracy']}")
    print(f"Edit distance vs referenca:    {s['edit_distance_reference']}")
    print(f"Saktesia vs referenca:         {s['accuracy_reference']}")
    print(f"Mutacione te futura:           {s['mutations_inserted']}")
    print(f"Mutacione te zbuluara sakte:   {s['mutations_correct']}")
    print(f"Alarme te rreme / te humbura:  {s['false_alarms']} / {s['missed']}")
    if not s["reference_available"]:
        print("Variantet:                     nuk llogariten pa reference")
    print(f"Koha e ekzekutimit:            {s['time_s']} s")
    print(f"Memoria maksimale (Python):    {s['peak_memory_mb']} MB")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="DNA genome reconstruction")
    ap.add_argument("--genome-len", type=int, default=333333)
    ap.add_argument("--reads", type=int, default=100000)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--k", type=int, default=12)
    ap.add_argument("--min-overlap", type=int, default=30)
    ap.add_argument("--max-error", type=float, default=0.08)
    ap.add_argument("--reads-file", help="Reads in FASTA/FASTQ, optionally .gz")
    ap.add_argument("--reference", help="Optional reference FASTA for variant reporting")
    args = ap.parse_args()

    corrected, final_rows, stats, output_contigs = run_pipeline(
        args.genome_len, args.reads, args.read_len, args.error_rate,
        args.k, args.min_overlap, args.max_error, args.reads_file, args.reference)
    print_report(stats)
    save_outputs(corrected, final_rows, output_contigs)
    print("U ruajten: reconstructed.fasta dhe variants_report.txt")
