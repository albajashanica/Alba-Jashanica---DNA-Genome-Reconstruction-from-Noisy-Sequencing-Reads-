# Rindërtimi i genomit të ADN-së nga reads të zhurmshme

## Për projektin

Ky projekt rindërton një sekuencë ADN-je nga reads të shkurtra, të parenditura dhe me gabime. Ai demonstron përdorimin e indeksimit me k-mer, zbulimin e overlap-eve, grafit të overlap-eve, montimit greedy dhe consensus-it. Projekti gjithashtu simulon variante dhe mat kohën e memorien për madhësi të ndryshme të dataset-it.

Të dhënat krijohen nga simulatori i projektit, prandaj rezultatet shërbejnë për demonstrim dhe vlerësim të algoritmit.

## Si punon

1. Krijohet një genom referencë dhe gjenerohen reads prej tij. Simulatori mund të shtojë gabime dhe të kthejë disa reads në orientimin reverse-complement.
2. Krijohen të dy orientimet për secilin read dhe ndërtohet indeks k-mer për të gjetur shpejt reads që mund të mbivendosen.
3. Kandidatët kontrollohen me edit distance për të vlerësuar cilësinë dhe gjatësinë e overlap-it.
4. Overlap-et paraqiten si lidhje të peshuara në graf. Montimi greedy i bashkon reads në contig-e.
5. Reads rreshtohen mbi contig-un kryesor. Consensus-i zgjedh bazat me mbështetjen më të madhe për të prodhuar sekuencën e korrigjuar.
6. Sekuenca krahasohet me referencën e simuluar për të gjetur dhe vlerësuar zëvendësime, inserime dhe fshirje.

## Skedarët

| Skedari | Përshkrimi |
|---|---|
| `main.py` | Nis pipeline-in kryesor dhe paraqet statistikat. |
| `simulator.py` | Gjeneron genom, reads, orientime dhe gabime të simuluara. |
| `indexing.py` | Krijon orientimet e reads dhe indeksin k-mer. |
| `overlap.py` | Gjen kandidatë dhe verifikon overlap-et. |
| `graph.py` | Ndërton dhe analizon grafikun e overlap-eve. |
| `assembler.py` | Ndërton grafikun për montim dhe krijon contig-e me metodën greedy. |
| `consensus.py` | Rreshton reads, ndërton pileup dhe llogarit consensus-in. |
| `mutations.py` | Simulon mutacione dhe krahason rindërtimin me referencën. |
| `filtering.py` | Filtron variantet sipas mbështetjes nga reads. |
| `scalability.py` | Teston projektin me numra të ndryshëm reads dhe ruan matjet në CSV. |
| `plot_results.py` | Krijon grafikët nga `scalability_results.csv`. |
| `reconstructed.fasta` | Shembull i sekuencës së rindërtuar. |
| `variants_report.txt` | Shembull i raportit të varianteve. |
| `scalability_results.csv` | Matje të ruajtura nga testi i shkallëzueshmërisë. |
| `scalability.png` | Grafik i rezultateve të shkallëzueshmërisë. |

## Ekzekutimi

Kërkohet Python 3. Nga dosja e projektit, ekzekuto pipeline-in kryesor:

```bash
python main.py
```

Për të ndryshuar parametrat, për shembull:

```bash
python main.py --genome-len 2000 --reads 300 --read-len 100 --error-rate 0.01 --k 12 --min-overlap 30 --max-error 0.08
```

Programi shtyp raportin në terminal dhe ruan sekuencën te `reconstructed.fasta`, ndërsa variantet i ruan te `variants_report.txt`.

## Testimi i shkallëzueshmërisë

Për të matur performancën me disa madhësi dataset-i:

```bash
python -u scalability.py 5000 10000 20000 --out scalability_results.csv
```

Kjo komandë i ekzekuton testet me 5,000, 10,000 dhe 20,000 reads dhe i ruan rezultatet te `scalability_results.csv`. Nëse ky fajll ekziston, rezultatet e mëparshme mbishkruhen. Për një test më të shpejtë që mat vetëm fazat e indeksimit dhe montimit, shto `--no-consensus`:

```bash
python -u scalability.py 5000 10000 20000 --no-consensus --out scalability_results.csv
```

CSV-ja përfshin kohën, memorien, numrin e kandidatëve dhe overlap-eve, si dhe përqindjen e genomit të mbuluar nga contig-u më i gjatë. Kur përdoret `--no-consensus`, statistikat që varen nga consensus-i nuk llogariten.

Për të krijuar grafikun nga `scalability_results.csv`:

```bash
python plot_results.py
```

Ky skript lexon pikërisht `scalability_results.csv` dhe ruan grafikun te `scalability.png`. Për të grafikuar rezultate nga një CSV tjetër, vendose atë me emrin `scalability_results.csv` para ekzekutimit.

## Rezultatet e raportuara

Pipeline-i raporton numrin e reads të përdorura, kandidatët dhe overlap-et e gjetura, cilësinë mesatare të overlap-eve, gjatësinë e rindërtimit, edit distance-in, saktësinë kundrejt mostrës, variantet e zbuluara, kohën e ekzekutimit dhe memorien maksimale të matur nga Python-i.

## Kufizime

Montimi greedy është heuristik dhe mund të mos zgjidhë saktë rajonet e përsëritura ose overlap-et e paqarta. Rezultatet varen nga parametrat dhe nga të dhënat e simuluara. Pipeline-i i tanishëm demonstron metodat mbi të dhëna të gjeneruara nga projekti; nuk është mjet i validuar për analiza klinike.



