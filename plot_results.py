import csv
import matplotlib.pyplot as plt

from matplotlib.ticker import NullFormatter

with open("scalability_results.csv") as f:
    rows = list(csv.DictReader(f))

reads = [int(r["reads"]) for r in rows]
naive = [int(r["naive"]) for r in rows]
cands = [int(r["candidates"]) for r in rows]
time_s = [float(r["t_total"]) for r in rows]
mem = [float(r["memory_mb"]) for r in rows]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

axes[0].loglog(reads, naive, "o-", label="Krahasime naive (n^2)")
axes[0].loglog(reads, cands, "s-", label="Kandidate me k-mer index")
axes[0].set_title("Numri i krahasimeve")
axes[0].set_xlabel("Numri i reads")
axes[0].set_ylabel("Krahasime (shkalle logaritmike)")
axes[0].legend()
axes[0].grid(True, which="both", alpha=0.3)

ideal_t = [time_s[0] * r / reads[0] for r in reads]
axes[1].plot(reads, time_s, "o-", label="Koha e matur")
axes[1].plot(reads, ideal_t, "--", label="Rritje lineare ideale")
axes[1].set_title("Koha e ekzekutimit")
axes[1].set_xlabel("Numri i reads")
axes[1].set_ylabel("Sekonda")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

ideal_m = [mem[0] * r / reads[0] for r in reads]
axes[2].plot(reads, mem, "o-", label="Memoria e matur")
axes[2].plot(reads, ideal_m, "--", label="Rritje lineare ideale")
axes[2].set_title("Memoria maksimale")
axes[2].set_xlabel("Numri i reads")
axes[2].set_ylabel("MB")
axes[2].legend()
axes[2].grid(True, alpha=0.3)

axes[0].set_xticks(reads)
axes[0].set_xticklabels(reads)
axes[0].xaxis.set_minor_formatter(NullFormatter())

plt.tight_layout()
plt.savefig("scalability.png", dpi=150)
print("U ruajt: scalability.png")
plt.show()