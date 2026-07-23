jest zwalidowany, ręcznie przetestowany pipeline (katalog test/) na jednym przykładzie z mimb_chapter/knowledge_bias/outputs: prepare_bias.py nadpisuje bias do map AD4, a vina --scoring ad4 --maps <prefix> doku­je z użyciem tych map bez potrzeby autodock4/DPF.

to działa, ale przykład miał już gotowe, wcześniej policzone mapy dołączone do repo

teraz celem jest benchmark metodologicznie kompletny, 3 nowe struktury krystalograficzne pobrane z PDB (benchmark/bialka/{6JQR,3CS9,5N9R}/*.cif + ligand *.sdf), dla których trzeba przygotować cały process od zera receptor, ligand, mapy AD4, miejsce biasu

pozniej porównać docking konwencjonalny vs biasowany w sposób

# WAŻNE
jak rozpoznać, które wody są potrzebne strukturalnie? - #TODO
    - cała wode usunąć(algorytm konkretny)
    jakie elementy chcemy zostawić do dokowania? - #TODO
    - NIC(zostawiamy tylko jony)
czy moge w kazdym miec jedna jednostke pomimo, ze pliki maja rozne jednostki? - #TODO
- jedno
czym jest luka w PDB, czy to jest luka w kieszeni? 

zainstalowanie autogrid4, #zrobione




prepare_receptor4.py i prepare_gpf4.py (Python2/AutoDockTools) są już dostępne w /opt/chimera/lib/python2.7/site-packages/AutoDockTools/Utilities24/ — używają dokładnie tego samego, uznanego protokołu AD4 co reszta pipeline'u.

- jedna wersja przez prepare_receptor4 i prepare_gpf4
    to zostało normalnie zastąpione przez Meeko od Forli Lab
    
- jedna przez skrypty Meeko + RDkit
    przygotowanie liganda (SDF bez H do PDBQT z typami AD4). Meeko + RDKit zamiast archaicznego prepare_ligand4.py w pythonie 3. RMSD symetryczny (potrzebny do oceny jakości pozy) policzymy przez RDKit, po odtworzeniu porządku wiązań z oryginalnym SDF. - #TODO
    naiwny RMSD po kolejności atomów byłby metodologicznie błędny dla symetrycznych podstawników (np. pierścienie).


benchmark typu redocking / pose-recovery. miejsce biasu dla każdego celu wyprowadzamy z własnych, znanych oddziaływań ligand–receptor tej samej struktury krystalicznej (bo mamy tylko jeden ligand na cel). 


bias / conv

PYTANIA DO HIPOTEZY
1. "czy podpowiedzenie znanej interakcji poprawia odzyskanie poprawnej pozy", żeby to nie była tylko korelacja a mechanizm 
2. "czy dodanie biasu gdziekolwiek w miejscu kieszeni a nie konkretnego wiązania nie powoduje jakieś zbieżności" dodatkowo
3. "czy to miejsce jest prawdopodobne ale błędne, nie powoduje problemów"

ZAŁOŻENIA:

1. Identycznie dla 3 celów (przygotowanie receptora, liganda, definicja pudełka), jedyna różnica między "conventional" i "biased" to obecność biasu.

2. Bias site wyznaczany w jaki sposób? #TODO 
    - programowo: geometryczne kryterium kontaktu donor/akceptor ligand a receptor struktury referencyjnej (odległość heteroatom–heteroatom < 5 Å + do sprawdzenia "rozsądny kąt"), udokumentowane dla każdego celu (który atom, która reszta, jaka odległość);
    na oko: ... #TODO; X

3. Ten sam Vset/r dla wszystkich 3 celów, trzeba sprawdzić konkretne cele #TODO, - 13 kcal/mol + 1.5 = - 11.5 kcal/mol. wyjasnienie z jednostka energi, ktora jest najmniejsza jako wiazanie wodorowe
    - normalnie (-1.5 kcal/mol, r=1.0 Å, jak w istniejącym przykładzie) - (acc, don, map, <aro>) - kilka jednocześnie, czy tylko jeden na przykład i tyle przykładów ile możliwych atomów, ideal interactions sites(pliku, testowanie kazdego), jak to jest zrobione
    - wymyślone przez nas( ... ) X

4. te same parametry(receptor(ligand +- r(15A), (r(10A)), r(5A), r(3A), ligand, box center, box size(dwie metody do sprawdzenia), exhaustiveness, num_modes, energy_range, seed, wersja Vina, liczba CPU).
    - raczej chcemy prostopadloscian a nie szescian. - do poczytania


5. Wspólna, sparowana pula seedów(6) dla conventional i biased tego samego celu (ten sam seed w obu warunkach), umożliwia to test parowany (np. Wilcoxon signed-rank), silniejszy statystycznie niż niesparowany. 
    - ile powtórzeń na test, z N=20 niezależnych powtórzeń na cel na warunek (2 warunki x 3 cele x 20 = 120 przebiegów Vina). 
        - co jest dokładnie parą statystyczną w tym przypadku, żeby to miało znaczenie ( .. ) do ustalenia
    - Metryka pierwszorzędna, symetryczny RMSD do natywnej pozy krystalograficznej, próg sukcesu ≤ 2.0 Å (konwencja redocking z literatury CASF/PDBbind) ALE #TODO(CZY TO PRAWDA I WARTO TO SPRAWDZAĆ ANIŻELI WSZYSTKIE ŚREDNIE) raportowany dla top-1 i best-of-N.
    - dlaczego ma być próg sukcesu <= X A, dlaczego X ma być 2 
    - inne testy(RMSD, docking score, czas, #TODO ( ... )) - Wynik energetyczny (score) NIE jest metryką porównawczą między conventional i biased — bias sztucznie zmienia energię w jednym miejscu siatki, więc bezwzględne wyniki nie są współmierne (już to ustaliliśmy przy pierwszym teście). Score służy tylko jako diagnostyka, nie jako "dowód" lepszego dockingu.

6. ilość seedów to jest do ustalenia, żeby znaleźć kilka pozycji i potem zrobić uśrednienie jako wynik wszystkiego oraz najlepsza pozycja z tych n seedów - #TODO jak zrobić indukcje z innych eksperymentów, zeby to miało znaczenie statystyczne

7. różne algorytmy przeszukiwania, AD4 i Vina, #TODO(X), trzeba przetestować oba silniki testowania, czy coś lepiej działa z Biasem
8. wszystkie kroki spięte w jeden config, skrypt i wersje narzędzi, żeby to było powtarzalne

9. kontrola negatywna? #TODO ( ... )



- metrka np. dipol ligandu a dopol kieszeni
- do konkretnego aminokwasu dokowanie(w molsmithie)

Nie w tym zakresie: kontrola negatywna z biasem w złym/losowym miejscu (test specyficzności) - rozszerzenie po tym, jak podstawowy benchmark zadziała, bo podwaja liczbę przebiegów i wymaga dodatkowej logiki wyboru miejsca pułapki.

  config.yaml                        # cele, seedy, padding pudełka, Vset/r, progi
  code/
    prepare_receptor.py           # CIF -> PDB (strip HETATM poza białkiem) -> prepare_receptor4.py -> receptor.pdbqt
    prepare_ligand.py             # SDF -> RDKit AddHs (pH~7.4) -> Meeko -> ligand.pdbqt (+ referencyjny mol do RMSD)
    extract_bias_site.py          # natywny ligand + receptor -> kontakty H-bond wg kryterium geometrycznego -> bias.bpf + raport tekstowy wyboru
    make_grid.py                  # prepare_gpf4.py -> receptor.gpf -> autogrid4 -> mapy .map
    run_docking.py                # pętla: cel x {conventional,biased} x seed -> prepare_bias.py (-b -g) + vina --scoring ad4
    score_rmsd.py                 # docked poses vs natywny ligand -> CalcRMS, success@2A, CSV per-run, statystyki, csv, wykres
  results/<PDBID>/{conventional,biased}/seed_<N>/ (out.pdbqt, vina.log)
  report/results.csv plots/