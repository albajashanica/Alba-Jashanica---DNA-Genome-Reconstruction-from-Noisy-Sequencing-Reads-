# Rindërtimi i genomit të ADN-së nga reads të zhurmshme

## Për projektin

Ky projekt rindërton një sekuencë ADN-je nga reads të shkurtra, të parenditura dhe me gabime. Ai demonstron përdorimin e indeksimit me k-mer, zbulimin e overlap-eve, grafit të overlap-eve, montimit greedy dhe consensus-it. Projekti gjithashtu simulon variante dhe mat kohën e memorien për madhësi të ndryshme të dataset-it.

Projekti mbështet si të dhëna të simuluara për vlerësim, ashtu edhe reads nga skedarë FASTA/FASTQ (përfshirë `.gz`). Një referencë FASTA është opsionale; pa të, raportohen rezultatet e montimit, por nuk bëhet thirrje për variante.

## Si punon

1. Reads merren nga një skedar FASTA/FASTQ ose gjenerohen nga simulatori. Simulatori mund të shtojë gabime dhe orientime reverse-complement.
2. Krijohen të dy orientimet për secilin read dhe ndërtohet indeks k-mer për të gjetur shpejt reads që mund të mbivendosen.
3. Kandidatët kontrollohen me edit distance për të vlerësuar cilësinë dhe gjatësinë e overlap-it.
4. Overlap-et paraqiten si lidhje të peshuara në graf. Montimi greedy zgjedh lidhjet më të forta pa përdorur një read dy herë dhe shmang ciklet.
5. Reads rreshtohen mbi contig-un kryesor. Consensus-i zgjedh bazat me mbështetjen më të madhe për të prodhuar sekuencën e korrigjuar.
6. Kur jepet referencë, sekuenca krahasohet me të për të gjetur dhe vlerësuar zëvendësime, inserime dhe fshirje.

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
python main.py --genome-len 333333 --reads 100000 --read-len 100 --error-rate 0.01 --k 12 --min-overlap 30 --max-error 0.08
```

Këta janë edhe parametrat e parazgjedhur: `python main.py` nis me 100,000 reads.

Për të montuar reads reale nga FASTA ose FASTQ. Kur jepet FASTQ, quality scores peshojnë votat e bazave në consensus:

```bash
python main.py --reads-file reads.fastq --k 15 --min-overlap 40 --max-error 0.08
```

Për të raportuar variacionet kundrejt një reference FASTA:

```bash
python main.py --reads-file reads.fastq --reference reference.fasta --k 15 --min-overlap 40 --max-error 0.08
```

Reads duhet të përmbajnë vetëm bazat A, C, G dhe T. Kur ka shumë reads ose referenca është e gjatë, zgjidh `k` dhe pragjet sipas gjatësisë së leximeve dhe mbulimit.

Programi shtyp raportin në terminal dhe ruan të gjitha contig-et në `reconstructed.fasta`; contig-u i parë është consensus-i i korrigjuar. Variantet ruhen te `variants_report.txt` dhe krahasohen në contig-un kryesor.

## Testimi i shkallëzueshmërisë

Për të matur performancën me disa madhësi dataset-i:

```bash
python -u scalability.py 25000 50000 100000 --no-consensus --out scalability_results.csv
```

Kjo komandë mat indeksimin, gjenerimin e kandidatëve, overlap-et dhe montimin për 25,000, 50,000 dhe 100,000 reads. `--no-consensus` kapërcen consensus-in dhe zbulimin e varianteve në këtë test të madh. Rezultatet e fundit të matura ruhen te `scalability_results.csv`; madhësia 100,000 mund të marrë shumë minuta dhe kërkon memorie të konsiderueshme. Pa `--no-consensus`, llogariten edhe consensus-i dhe variantet, por vetëm për madhësi të vogla.

Për një test më të vogël që përfshin consensus-in dhe zbulimin e varianteve:

```bash
python -u scalability.py 300 600 1200 --out scalability_full.csv
```

CSV-ja përfshin kohën, memorien, numrin e kandidatëve dhe overlap-eve, si dhe përqindjen e genomit të mbuluar nga contig-u më i gjatë. Për 10,000 reads e më pak memoria matet me `tracemalloc`; për dataset-et më të mëdha matet memoria kulmore e procesit, për të shmangur overhead-in e gjurmimit të çdo alokimi. Indeksi heq k-mer-at me mbi 200 shfaqje për të kufizuar listat e kandidatëve; kjo mund të ulë ndjeshmërinë në rajone shumë të përsëritura. Kur përdoret `--no-consensus`, statistikat që varen nga consensus-i nuk llogariten.

Matjet e ekzekutuara në këtë projekt:

| Reads | Kandidatë | Overlap-e | Contig-u më i gjatë | Koha (s) | Memoria kulmore (MB) |
|---:|---:|---:|---:|---:|---:|
| 25,000 | 1,302,892 | 1,046,215 | 56.0% | 368.03 | 744.7 |
| 50,000 | 2,697,710 | 2,092,990 | 79.5% | 541.17 | 1,287.6 |
| 100,000 | 5,553,648 | 4,177,526 | 99.8% | 1,031.69 | 2,138.1 |

Këto pika përdorin `--no-consensus`; koha dhe memoria nuk përfshijnë consensus-in dhe thirrjen e varianteve.

Për të krijuar grafikun nga `scalability_results.csv`:

```bash
python plot_results.py
```

Ky skript lexon pikërisht `scalability_results.csv` dhe ruan grafikun te `scalability.png`. Për të grafikuar rezultate nga një CSV tjetër, vendose atë me emrin `scalability_results.csv` para ekzekutimit.

## Rezultatet e raportuara

Pipeline-i raporton numrin e reads të përdorura, kandidatët dhe overlap-et e gjetura, cilësinë mesatare të overlap-eve, gjatësinë e rindërtimit, edit distance-in, saktësinë kundrejt mostrës, variantet e zbuluara, kohën e ekzekutimit dhe memorien maksimale të matur nga Python-i.

## Kufizime

Montimi greedy është heuristik dhe mund të mos zgjidhë saktë rajonet e përsëritura ose overlap-et e paqarta. Matjet e ruajtura deri në 100,000 reads përdorin `--no-consensus`, prandaj ato matin indeksimin, overlap-et dhe montimin, jo consensus-in apo zbulimin e varianteve në atë madhësi. K-mer-at që shfaqen mbi 200 herë hiqen nga indeksi, gjë që mund të ulë ndjeshmërinë në rajone të përsëritura. Saktësia dhe zbulimi i varianteve kërkojnë referencë të njohur; rezultatet varen nga parametrat dhe mbulimi. Projekti është për demonstrim dhe nuk është mjet i validuar për analiza klinike.



