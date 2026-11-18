import csv
import matplotlib.pyplot as plt

with open("scalability_results.csv") as f:
    rows = list(csv.DictReader(f))

reads = [int(r["reads"]) for r in rows]
naive = [int(r["naive"]) for r in rows]
cands = [int(r["candidates"]) for r in rows]
time_s = [float(r["t_total"]) for r in rows]
mem = [float(r["memory_mb"]) for r in rows]