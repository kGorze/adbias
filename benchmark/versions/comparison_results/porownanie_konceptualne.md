### Równoważność konceptualna new względem legacy:
- Cele: 6JQR, 3CS9, 5N9R.
- Warunki: `conventional` i `biased`.
- Wspólne seedy: 1–5.
- Liczba uczciwie sparowanych przebiegów: 30 na wersję, 60 łącznie.
- Próg sukcesu: RMSD ≤ 2 Å.
- 5 seedow, 5 seedow 

Klasyfikacja sukces/porażka jest zgodna w 12/12 porównań (3 cele x 2 warunki x top-1/best-of-N).

| cel | warunek | legacy top-1/best success | new top-1/best success | wniosek |
|---|---|---:|---:|---|
| 6JQR | conventional | 0% / 0% | 0% / 0% | zgodny: porażka |
| 6JQR | biased | 0% / 0% | 0% / 0% | zgodny: porażka |
| 3CS9 | conventional | 100% / 100% | 100% / 100% | zgodny: sukces |
| 3CS9 | biased | 100% / 100% | 100% / 100% | zgodny: sukces |
| 5N9R | conventional | 100% / 100% | 100% / 100% | zgodny: sukces |
| 5N9R | biased | 100% / 100% | 100% / 100% | zgodny: sukces |

Kontrola metryki na tych samych 120 pozach top-1 z `new` daje medianę i maksimum `|RDKit − obrms|` równe 0.000 Å. Różnice między wersjami nie wynikają więc ze zmiany kalkulatora RMSD w tym zbiorze.

Logika wyboru biasu daje identyczne miejsce, typ oraz odległość w 3/3 celach:
| cel | kontakt wybrany przez obie wersje | typ | odległość |
|---|---|---|---:|
| 6JQR | ligand O (OA) ↔ CYS694/N (N) | acc | 2.80 Å |
| 3CS9 | ligand N (N) ↔ THR315/OG1 (OA) | don | 2.92 Å |
| 5N9R | ligand O (OA) ↔ TYR465/OH (OA) | acc | 2.47 Å |

Pozostałe parametry eksperymentu — cele, warunki, seedy, liczba mód, energy range, spacing, padding, głębokość i promień biasu, cutoff kontaktu oraz próg sukcesu RMSD — są zgodne dla danych.

W `new` Meeko usuwa część niekompletnych reszt w 6JQR i 3CS9.

### Jawna rozbieżność i granica argumentacji:

Kierunek niewielkiej zmiany mediany top-1 po włączeniu biasu nie zgadza się dla 3CS9:

| cel | legacy: biased − conventional | new: biased − conventional | zgodność kierunku |
|---|---:|---:|---|
| 6JQR | −0.361 Å | −0.597 Å | tak, poprawa |
| 3CS9 | −0.009 Å | +0.048 Å | nie |
| 5N9R | −0.008 Å | −0.021 Å | tak, poprawa |

- `comparison_runs.csv` — wyniki per wersja/cel/warunek/seed; wiersze wspólne identyfikuje ta sama trójka `pdbid`, `condition`, `seed`.
- `comparison_summary.csv` — agregaty wyłącznie dla wspólnych seedów.

- `../compare.py` — skrypt generujący oba CSV-y i kontrolę RDKit vs `obrms`.
- `../new/results/<PDBID>/receptor_report.txt` — kontrola usuniętych reszt i ich odległości od liganda.
- `../{legacy,new}/results/<PDBID>/bias_report.txt` — źródłowe raporty wyboru miejsca biasu.

### Odtworzenie:

```bash
/home/kgorzelanczyk/miniforge3/envs/vs/bin/python benchmark/versions/compare.py
```
