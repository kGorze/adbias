jest zwalidowany, ręcznie przetestowany pipeline (katalog test/) na jednym przykładzie z mimb_chapter/knowledge_bias/outputs: prepare_bias.py nadpisuje bias do map AD4, a vina --scoring ad4 --maps <prefix> doku­je z użyciem tych map bez potrzeby autodock4/DPF.

to działa, ale przykład miał już gotowe, wcześniej policzone mapy dołączone do repo

teraz celem jest benchmark metodologicznie kompletny, 3 nowe struktury krystalograficzne pobrane z PDB (benchmark/bialka/{6JQR,3CS9,5N9R}/*.cif + ligand *.sdf), dla których trzeba przygotować cały process od zera receptor, ligand, mapy AD4, miejsce biasu

pozniej porównać docking konwencjonalny vs biasowany w sposób

jak rozpoznać, które wody są potrzebne strukturalnie? - #TODO



Braki narzędziowe zidentyfikowane i rozwiązane:

autogrid4 nie jest zainstalowany, ale jest dostępny: conda install -c conda-forge autogrid=4.2.9. Bez tego nie da się policzyć map AD4 dla nowych struktur (Vina w trybie ad4 tylko czyta gotowe mapy, nigdy ich nie liczy).
prepare_receptor4.py i prepare_gpf4.py (Python2/AutoDockTools) są już dostępne w /opt/chimera/lib/python2.7/site-packages/AutoDockTools/Utilities24/ — używają dokładnie tego samego, uznanego protokołu AD4 co reszta pipeline'u.
Do przygotowania liganda (SDF bez H → PDBQT z typami AD4) użyjemy Meeko + RDKit (env bo-mmpbsa-docking, ma obie biblioteki) zamiast archaicznego prepare_ligand4.py — nowocześniejsze, mniej kruche, Python3.
RMSD symetryczny (potrzebny do oceny jakości pozy) policzymy przez rdkit.Chem.rdMolAlign.GetBestRMS po odtworzeniu porządku wiązań z AssignBondOrdersFromTemplate (referencją jest oryginalny SDF) — naiwny RMSD po kolejności atomów byłby metodologicznie błędny dla symetrycznych podstawników (np. pierścienie).
Zakres i uczciwe ramy metodologiczne (kluczowe dla "bez zarzutu")
To jest benchmark typu redocking / pose-recovery, nie ślepy virtual screening: miejsce biasu dla każdego celu wyprowadzamy z własnych, znanych oddziaływań ligand–receptor tej samej struktury krystalicznej (bo mamy tylko jeden ligand na cel — nie ma z czego zrobić "ślepego" testu). To jest dokładnie ten sam projekt eksperymentu co w oryginalnej publikacji AutoDock Bias — pytanie brzmi "czy podpowiedzenie znanej interakcji poprawia odzyskanie poprawnej pozy", a NIE "czy metoda przewiduje nieznane wiązanie". Ta rama musi być jawnie opisana w raporcie końcowym, żeby uniknąć zarzutu obiegu (circularity).

Zasady rygoru:

Identyczny protokół dla 3 celów (przygotowanie receptora, liganda, definicja pudełka) — jedyna różnica między "conventional" i "biased" to obecność biasu.
Bias site wyznaczany programowo, nie na oko: geometryczne kryterium kontaktu donor/akceptor ligand↔receptor (odległość heteroatom–heteroatom ≤ 3.5 Å + rozsądny kąt), udokumentowane dla każdego celu (który atom, która reszta, jaka odległość).
Ten sam Vset/r dla wszystkich 3 celów (-1.5 kcal/mol, r=1.0 Å — jak w istniejącym przykładzie) — inaczej można zarzucić dobieranie parametrów pod wynik.
Wspólna, sparowana pula seedów dla conventional i biased tego samego celu (ten sam seed w obu warunkach) → pozwala na test parowany (Wilcoxon signed-rank), silniejszy statystycznie niż niesparowany.
N=20 niezależnych powtórzeń na cel na warunek (2 warunki × 3 cele × 20 = 120 przebiegów Vina — z tym ligandem to rzędu minut, nie godzin).
Metryka pierwszorzędna: symetryczny RMSD do natywnej pozy krystalograficznej, próg sukcesu ≤ 2.0 Å (konwencja redocking z literatury CASF/PDBbind), raportowany dla top-1 i best-of-N.
Wynik energetyczny (score) NIE jest metryką porównawczą między conventional i biased — bias sztucznie zmienia energię w jednym miejscu siatki, więc bezwzględne wyniki nie są współmierne (już to ustaliliśmy przy pierwszym teście). Score służy tylko jako diagnostyka, nie jako "dowód" lepszego dockingu.
Pre-flight check luk w sekwencji (np. 6JQR 29 niemodelowanych reszt) względem odległości do centroidu liganda — jeśli luka jest blisko kieszeni, flagować to jawnie w raporcie jako ograniczenie.
Wszystkie kroki spięte jednym skryptem/configiem, wersje narzędzi (vina --version, autogrid4 --version, meeko/rdkit) zapisane do logu — powtarzalność.
Nie w tym zakresie (świadomie odłożone): kontrola negatywna z biasem w złym/losowym miejscu (test specyficzności) — to naturalne rozszerzenie po tym, jak podstawowy benchmark zadziała, bo podwaja liczbę przebiegów i wymaga dodatkowej logiki wyboru miejsca "decoy". Zaproponuję to jako Fazę 2 po pierwszych wynikach.

Struktura plików do stworzenia
benchmark/
  config.yaml                        # cele, seedy, padding pudełka, Vset/r, progi
  scripts/
    01_prepare_receptor.py           # CIF -> PDB (strip HETATM poza białkiem) -> prepare_receptor4.py -> receptor.pdbqt
    02_prepare_ligand.py             # SDF -> RDKit AddHs (pH~7.4) -> Meeko -> ligand.pdbqt (+ referencyjny mol do RMSD)
    03_extract_bias_site.py          # natywny ligand + receptor -> kontakty H-bond wg kryterium geometrycznego -> bias.bpf + raport tekstowy wyboru
    04_make_grid.py                  # prepare_gpf4.py -> receptor.gpf -> autogrid4 -> mapy .map
    05_run_docking.py                # pętla: cel x {conventional,biased} x seed -> prepare_bias.py (-b -g) + vina --scoring ad4
    06_score_rmsd.py                 # docked poses vs natywny ligand -> GetBestRMS, success@2A, CSV per-run
    07_aggregate_report.py           # statystyki (Wilcoxon sparowany), tabela zbiorcza, wykresy, summary.md
  results/<PDBID>/{conventional,biased}/seed_<N>/ (out.pdbqt, vina.log)
  report/summary.csv, summary.md, plots/
Każdy skrypt jest małym, pojedynczym krokiem CLI (wejście/wyjście = pliki), żeby dało się uruchomić i zweryfikować krok po kroku zanim poleci pełna pętla na 3×2×20 przebiegów.

Kolejność wykonania (plan działania)
conda install -c conda-forge autogrid=4.2.9 do dedykowanego env (np. nowy benchmark albo istniejący adbias) + zapisanie wersji do loga.
Napisanie i przetestowanie skryptów 01–04 na jednym celu — 3CS9 (najczystsza struktura: tylko HOH poza ligandem) — jako pilot, z ręczną inspekcją wygenerowanego receptor.pdbqt, receptor.gpf, map.
Pre-flight check luk sekwencyjnych (punkt 8 wyżej) dla wszystkich 3 celów — jeśli 6JQR ma problem, zdecydujemy wtedy (np. modeler loop albo udokumentowane ograniczenie).
Sanity check dockingu: pozycyjny "positive control" — dokowanie liganda startującego z natywnej pozy z --local_only (bez randomizacji) powinno dać RMSD ~0 Å; to potwierdza, że przygotowanie receptor/ligand/mapy jest poprawne, zanim zaufamy pełnym wynikom.
Rozszerzenie na 6JQR i 5N9R, uruchomienie skryptu 03 (bias site) i ręczna weryfikacja sensowności wybranego kontaktu dla każdego z 3 celów (krótki opis w raporcie: która reszta, jaki atom, jaka odległość).
Pełny przebieg 05 (120 dockingów), 06 (scoring), 07 (raport + statystyka).
Przegląd report/summary.md razem z użytkownikiem — dopiero potem ewentualna Faza 2 (kontrola negatywna).
Weryfikacja
Krok 4 powyżej (positive control RMSD≈0) jest twardym testem poprawności przygotowania każdego celu — brak tego = nie ufamy dalszym wynikom.
Po pełnym przebiegu: report/summary.md musi zawierać per-cel i zbiorczo — medianę/IQR RMSD, success rate @2Å (top-1 i best-of-N), wynik testu Wilcoxona, oraz jawny opis wybranego miejsca biasu (reszta/atom/odległość) dla każdego z 3 celów.
Logi surowe (results/.../vina.log) zachowane, żeby każdy wynik dało się prześledzić do konkretnego przebiegu i seeda.