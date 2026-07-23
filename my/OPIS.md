## Najważniejszy werdykt

Plan ma dobry rdzeń techniczny, ale obecnie miesza cztery różne pytania badawcze:

1. Czy implementacja C++ reprodukuje zachowanie AutoDock Bias?
2. Czy bias wyprowadzony z natywnej pozy liganda poprawia redocking?
3. Czy poprawa jest specyficzna dla prawidłowej interakcji, a nie dla dowolnego minimum energii w kieszeni?
4. Czy bias wyprowadzony bez znajomości pozy badanego liganda poprawia docking nowych ligandów?

Te pytania wymagają osobnych eksperymentów. Trzy struktury i 20 seedów wystarczą do pilota technicznego i mechanistycznego. Nie wystarczą do wykazania generalizacji metody na populację receptorów lub ligandów.

Warto też sprostować punkt odniesienia: katalog `mimb_chapter` jest tutorialem i przykładem regresyjnym, ale oryginalna walidacja AutoDock Bias była szersza. Autorzy wykonali cross-docking 10 ligandów AmpC, użyli biasów pochodzących z symulacji mixed-solvent MD i oceniali zarówno odzyskanie pozy, jak i virtual screening. Za poprawną uznawali pozę z ciężkoatomowym RMSD poniżej 2 Å. ([OUP Academic][1])

Najlepsza koncepcja projektu brzmi zatem:

> Nowa implementacja powinna najpierw dokładnie odtwarzać AutoDock Bias, następnie oddzielnie wykazać mechanizm działania, specyficzność prawidłowego biasu oraz generalizację biasów wyprowadzonych bez dostępu do badanego liganda.

---

# 1. Cztery poziomy benchmarku

## Poziom 0: zgodność implementacyjna

Cel: udowodnić, że C++ poprawnie zastępuje `prepare_bias.py`.

Nie wykonuj tutaj jeszcze benchmarku dockingu. Porównuj bezpośrednio wartości map:

[
\Delta(g)=V_{\mathrm{C++}}(g)-V_{\mathrm{Python2}}(g)
]

dla każdego punktu siatki i każdego typu mapy.

AutoDock Bias dodaje do oryginalnej mapy odwróconą funkcję Gaussa:

[
V_{\mathrm{biased}}(g)=V_{\mathrm{original}}(g)+
V_{\mathrm{set}}
\exp\left(-\frac{|g-c|^2}{r^2}\right)
]

gdzie (V_{\mathrm{set}}<0). ([OUP Academic][1])

Ważna korekta Twojego zapisu: jeżeli oryginalna wartość mapy w centrum wynosi `-13 kcal/mol`, a `Vset=-1.5 kcal/mol`, to nowa wartość wynosi:

[
-13 + (-1.5)=-14.5\ \mathrm{kcal/mol}
]

Nie `-11.5 kcal/mol`.

Dla `Vset=-1.5 kcal/mol` i `r=1 Å`:

* w centrum: dodatkowo `-1.500 kcal/mol`,
* w odległości `r`: `-0.552 kcal/mol`,
* w odległości `1.5r`: około `-0.158 kcal/mol`,
* w odległości `2r`: około `-0.027 kcal/mol`.

Parametr `r` nie jest twardym promieniem odcięcia. Jest promieniem zaniku typu (e^{-1}). Po zapisaniu funkcji jako klasycznej Gaussowskiej wariancji odpowiada (\sigma=r/\sqrt{2}).

Testy zgodności powinny obejmować:

* bias dokładnie w punkcie siatki,
* bias pomiędzy punktami siatki,
* bias przy granicy mapy,
* kilka nakładających się biasów,
* `don`, `acc`, `aro` i dowolną wskazaną mapę,
* `Vset=0`,
* dodatnią karę, jeżeli nowa implementacja ma ją obsługiwać,
* miejsce całkowicie poza mapą — powinien być jawny błąd, nie ciche obcięcie.

To jest test zgodności oprogramowania, nie test naukowy. Wynik dockingu jest zbyt pośredni, żeby wykryć subtelny błąd indeksowania, położenia początku siatki albo interpretacji `r`.

---

## Poziom 1: oracle redocking

Tutaj bias jest wyprowadzany z pozy krystalograficznej tego samego liganda.

To jest poprawny eksperyment dla hipotezy:

> Czy podanie algorytmowi prawidłowej informacji o konkretnej interakcji zwiększa prawdopodobieństwo odnalezienia lub prawidłowego uszeregowania pozy natywnej?

Nie jest to jednak test predykcyjny. Używasz informacji z odpowiedzi do utworzenia wejścia. Należy to nazywać:

* `oracle bias`,
* `native-derived bias`,
* albo `upper-bound experiment`.

Nie nazywałbym tego ogólnie „receptor-derived bias”, ponieważ w tym wariancie informacja pochodzi z kompleksu receptor–badany ligand.

## Poziom 2: specyficzność i odporność

Tu trzeba wykazać, że działa właśnie prawidłowa interakcja, a nie dowolna sztuczna studnia energii.

Minimalne warunki:

1. `no_bias`,
2. `correct_oracle_bias`,
3. `matched_decoy_bias`.

Dodatkowe, bardzo wartościowe ablac­je:

4. `wrong_feature_type` — prawidłowa pozycja, ale błędny typ, np. donor zamiast akceptora,
5. `shifted_1A`,
6. `shifted_2A`,
7. `shifted_4A`,
8. `sham_bias` z `Vset=0`.

Kontrola negatywna nie może zostać odłożona, jeżeli chcesz twierdzić, że efekt jest specyficzny lub mechanistyczny. Bez niej możesz powiedzieć jedynie:

> Dodanie prawidłowo umieszczonego minimum energii poprawiło wyniki w eksperymencie oracle.

Nie możesz wtedy powiedzieć:

> Poprawa wynikała z rozpoznania konkretnej interakcji receptor–ligand.

## Poziom 3: generalizacja

To jest właściwy benchmark nowej metody.

Bias dla badanego liganda musi pochodzić z informacji niezależnej od jego pozy, na przykład:

* z innych ligandów tego receptora,
* z homologicznego receptora,
* z mutagenezy,
* z mixed-solvent MD,
* z map hotspotów,
* z apo-struktury,
* z konserwowanych oddziaływań w rodzinie białek,
* z przewidywanych miejsc wody lub metalu.

Najlepszy układ to leave-one-ligand-out lub leave-one-chemotype-out:

1. Dla receptora dostępnych jest kilka kompleksów.
2. Jeden ligand jest całkowicie odkładany jako testowy.
3. Bias wyprowadza się wyłącznie z pozostałych ligandów lub danych receptorowych.
4. Testowy ligand jest dokowany do receptora, którego wybór również nie zależy od jego pozy.
5. Podział wykonuje się po szkielecie chemicznym, a nie losowo po ligandach.

Dopiero ten poziom pozwala stwierdzić, że metoda pomaga w przewidywaniu nowych pozy.

---

# 2. Poprawione hipotezy

Twoje hipotezy warto zapisać formalnie.

## H1 — efekt oracle

Przy identycznym przygotowaniu układu i identycznym budżecie obliczeniowym prawidłowy bias zwiększa prawdopodobieństwo uzyskania fizycznie poprawnej pozy top-1 o RMSD ≤ 2 Å względem dockingu bez biasu.

## H2 — specyficzność przestrzenna

Efekt prawidłowego biasu jest większy niż efekt biasu umieszczonego w dopasowanym, ale błędnym miejscu tej samej kieszeni.

## H3 — specyficzność chemiczna

Bias o prawidłowej lokalizacji, lecz błędnym typie interakcji daje mniejszą poprawę niż bias o prawidłowej lokalizacji i typie.

## H4 — odpowiedź na perturbację

Skuteczność maleje w kontrolowany sposób wraz z przesuwaniem miejsca biasu od prawidłowej lokalizacji.

## H5 — generalizacja

Bias wyprowadzony bez użycia badanego liganda zwiększa skuteczność na ligandach i chemotypach niewykorzystanych do wyznaczania biasu.

## H6 — interakcja z silnikiem

Efekt biasu zależy od kombinacji algorytmu przeszukiwania, bazowej funkcji scoringowej i reprezentacji biasu.

H1–H4 są pytaniami mechanistycznymi dotyczącymi algorytmu. H5 jest pytaniem predykcyjnym. Nie należy ich łączyć w jeden test.

---

# 3. Co dokładnie jest jednostką statystyczną

To jest najważniejsza poprawka statystyczna.

Dwadzieścia seedów dla jednego kompleksu nie oznacza 20 niezależnych eksperymentów biologicznych. Są to repliki techniczne tego samego problemu receptor–ligand.

Masz zatem:

* 3 niezależne kompleksy,
* po 20 stochastycznych realizacji każdego kompleksu.

Nie masz (n=60) niezależnych kompleksów.

Traktowanie wszystkich seedów ze wszystkich kompleksów jako jednej próby w teście Wilcoxona byłoby pseudoreplikacją.

## Definicja pary

Prawidłowa para statystyczna to:

> ten sam kompleks, ta sama wersja receptora, ten sam mikrostan liganda, ten sam box, ten sam silnik, ten sam budżet obliczeniowy i ten sam seed; różni się wyłącznie warunek biasu.

Wspólny seed jest użyteczny jako blok typu common random numbers. Nie oznacza jednak, że oba dockingi podążają tą samą trajektorią. Dokumentacja Vina wskazuje, że dokładna powtarzalność wymaga identyczności wszystkich wejść i parametrów; nawet mała zmiana wejścia może zachowywać się podobnie do zmiany seedu. Zmiana mapy przez bias właśnie taką zmianą jest. ([autodock-vina.readthedocs.io][2])

## Co można zrobić z trzema kompleksami

Dla każdego kompleksu osobno raportuj:

* udział sukcesów top-1,
* przedział ufności udziału,
* medianę i rozkład RMSD,
* różnicę udziałów sukcesu:
  [
  \Delta p_c=p_{c,\mathrm{correct}}-p_{c,\mathrm{no\ bias}},
  ]
* wyniki sparowane seed po seedzie,
* liczbę par:

  * porażka → sukces,
  * sukces → porażka.

Można użyć dokładnego testu McNemara wewnątrz kompleksu, ale będzie on opisywał powtarzalność algorytmiczną dla tego jednego układu, a nie generalizację na inne receptory.

Przy trzech kompleksach nawet wszystkie trzy efekty w tym samym kierunku dają w dwustronnym teście znaków (p=0.25). To pokazuje, dlaczego trzy struktury należy nazywać pilotem, a nie benchmarkiem potwierdzającym.

## Docelowy model dla większego zbioru

Dla wielu kompleksów użyłbym hierarchicznego modelu logistycznego:

[
\operatorname{logit}P(Y_{csk}=1)=
\beta_0+
\beta_1 I(k=\mathrm{correct})+
\beta_2 I(k=\mathrm{decoy})+
u_c+
v_{cs}
]

gdzie:

* (Y_{csk}=1), jeżeli top-1 ma RMSD ≤ 2 Å i przechodzi kontrolę fizycznej poprawności,
* (u_c) jest efektem losowym kompleksu,
* (v_{cs}) jest efektem bloku kompleks–seed,
* (k) oznacza warunek.

Dla wielu ligandów na jeden target dodałbym osobne efekty losowe targetu i kompleksu zagnieżdżonego w targecie.

Raportuj przede wszystkim:

* różnicę ryzyka,
* iloraz szans,
* medianę sparowanej zmiany RMSD,
* 95-procentowe przedziały ufności,
* wyniki per target.

Nie tylko wartości (p).

---

# 4. Ile seedów

Seed nie zwiększa liczby niezależnych receptorów. Zmniejsza jedynie niepewność Monte Carlo dla danego receptor–ligand.

Dla udziału sukcesów najgorszy standardowy błąd przy (S=20) wynosi:

[
\sqrt{\frac{0.5(1-0.5)}{20}}\approx 0.112
]

czyli przybliżony 95-procentowy margines to około ±0.22. To dość szeroko, ale wystarcza do pilota.

Praktyczna interpretacja:

* 6 seedów: smoke test,
* 10–20 seedów: pilot i wykrywanie dużych efektów,
* około 40 seedów: dokładniejsze oszacowanie sukcesu pojedynczego kompleksu,
* dalszy budżet lepiej przeznaczać na nowe niezależne kompleksy niż na setki seedów tego samego układu.

Jeżeli masz budżet zbliżony do planowanych 120 przebiegów na silnik, lepszy eksperyment to:

[
3\ \text{kompleksy}\times
3\ \text{warunki}\times
13\ \text{seedów}=117
]

Warunki:

* `no_bias`,
* `correct_bias`,
* `matched_decoy`.

To daje znacznie więcej informacji niż:

[
3\times2\times20=120
]

bez kontroli negatywnej.

Jeżeli produkcyjny algorytm ma działać jako „uruchom 20 seedów i wybierz najlepszy wynik”, to cały pakiet 20 seedów jest jedną realizacją algorytmu. Wtedy trzeba powtarzać niezależne pakiety, a nie traktować poszczególnych seedów jako wynik końcowy.

Dla pełnego benchmarku liczebność należy wyznaczyć na podstawie wariancji z pilota. Orientacyjnie, dla sparowanego efektu ciągłego o standaryzowanej wielkości (d_z=0.5) potrzeba około 34 niezależnych kompleksów przy mocy 80%; dla (d_z=0.3) około 90. Dla głównej metryki binarnej i modelu hierarchicznego właściwsza będzie symulacja mocy na podstawie wyników pilota.

---

# 5. Metryka pierwszorzędna

Proponuję:

[
Y=1[
\mathrm{RMSD}_{\mathrm{top1}}\leq 2.0\ \text{Å}
\land
\mathrm{pose_valid}
]
]

Czyli sukces wymaga jednocześnie:

1. symetrycznie poprawnego, ciężkoatomowego RMSD top-1 ≤ 2 Å,
2. braku poważnych błędów chemicznych lub sterycznych.

Próg 2 Å jest uznaną konwencją dockingową i był użyty także w oryginalnym AutoDock Bias, ale nie jest prawem fizycznym. Dlatego należy równolegle raportować ciągły RMSD i analizę czułości dla progów 1, 2 i 3 Å. ([OUP Academic][1])

RMSD sam w sobie nie wystarcza. PoseBusters pokazał, że natywopodobne RMSD może współistnieć z nieprawidłową stereochemią, geometrią wiązań, planarity lub kolizjami białko–ligand. ([PubMed][3])

## Krytyczna poprawka dotycząca RDKit

Do oceny dockingowej użyj `CalcRMS`, nie `GetBestRMS`.

`CalcRMS`:

* liczy RMSD „in place”,
* nie dopasowuje ponownie liganda do referencji,
* uwzględnia symetrię,
* jest w dokumentacji wskazany jako odpowiedni do porównywania pozy dockingowych z ligandem krystalograficznym. ([rdkit.org][4])

`GetBestRMS` wykonuje optymalne dopasowanie przestrzenne. W dockingu byłoby to metodologicznie niebezpieczne, ponieważ może usunąć translację i rotację błędnie umieszczonego liganda.

Procedura:

1. Zachowaj wspólny układ współrzędnych receptora.
2. Odtwórz ligand z prawidłowymi bond orders i stereochemią.
3. Usuń wodory z obliczenia RMSD.
4. Użyj jawnego mapowania atomów wynikającego z oryginalnego SDF.
5. Uwzględnij automorfizmy i grupy równoważne, np. karboksylany i symetryczne podstawniki.
6. Nie pozwalaj algorytmowi mapowania na zmianę tautomeru lub formalnego ładunku.

## Metryki drugorzędne

Raportowałbym:

* ciągły `top1_RMSD`,
* `success@1A`, `success@2A`, `success@3A`,
* `best-of-K RMSD`, przy stałym (K),
* pozycję pierwszej natywopodobnej pozy w rankingu,
* reciprocal rank pierwszej poprawnej pozy,
* spełnienie wskazanej interakcji,
* zgodność interaction fingerprint,
* kolizje białko–ligand,
* strain liganda,
* liczbę klastrów i różnorodność pozy,
* czas CPU i wall time,
* liczbę ewaluacji funkcji celu,
* zużycie pamięci.

`Best-of-K` nie powinno być metryką pierwszorzędną, jeżeli realne zastosowanie zwraca top-1. To metryka zdolności generowania pozy, nie zdolności prawidłowego rankingu.

---

# 6. Rozdzielenie search failure i scoring failure

Nowa implementacja może zrobić coś, czego stary AdBias nie eksponuje jasno.

Dla każdej pozy zapisuj trzy energie:

[
E_{\mathrm{search}}=E_{\mathrm{base}}+E_{\mathrm{bias}}
]

oraz oddzielnie:

* `base_score` — wynik na niezmodyfikowanych mapach,
* `bias_score` — wkład biasu,
* `search_total` — suma używana podczas przeszukiwania.

Dodatkowo wykonaj tani rescore wszystkich znalezionych pozy na czystych mapach.

Daje to trzy protokoły:

1. `base_search + base_rank`,
2. `biased_search + biased_rank`,
3. `biased_search + base_rank`.

Porównanie 1 z 3 mierzy, czy bias poprawia generowanie pozy.

Porównanie 2 z 3 mierzy, czy bias jest potrzebny także do końcowego rankingu.

To umożliwia klasyfikację porażek:

* **search failure** — w top-K nie wygenerowano poprawnej pozy,
* **ranking failure** — poprawna poza powstała, ale nie została wybrana,
* **constraint underspecification** — bias został spełniony, ale cała poza jest błędna,
* **chemistry failure** — poza jest chemicznie lub sterycznie nieprawidłowa,
* **preparation/reference failure** — problem z protonacją, mapowaniem atomów lub strukturą.

Masz rację, że bezwzględnego biased score nie wolno bezpośrednio porównywać z conventional score. To są różne funkcje celu. Porównywalny jest dopiero `base_score` przeliczony dla wszystkich pozy na tej samej, niezmodyfikowanej funkcji.

---

# 7. Wyznaczanie miejsca biasu

Kryterium „heteroatom–heteroatom < 5 Å” jest za luźne jako definicja wiązania wodorowego. Może służyć jedynie jako wstępne wyszukanie kandydatów.

Dla oracle-pilota proponuję następującą procedurę:

1. Przygotować prawidłową protonację receptora i liganda.
2. Wyszukać wszystkie chemicznie zgodne pary donor–akceptor.
3. Jako prerejestrowany punkt startowy użyć:

   * donor–akceptor ≤ 3.5 Å,
   * wodór–akceptor ≤ około 2.6 Å,
   * kąt D–H···A ≥ 120°,
4. Odrzucić atomy z niską occupancy, alternatywną konformacją lub słabą gęstością.
5. Uszeregować kontakty według jakości geometrii, kompletności atomów i jednoznaczności chemicznej.
6. Wybrać jedno miejsce przed uruchomieniem dockingu.
7. Zapisać pełny raport wyboru.

PLIP wykrywa regułowo m.in. wiązania wodorowe, oddziaływania hydrofobowe, stacking, salt bridges, water bridges i halogen bonds. Arpeggio rozszerza to o szerszy zestaw kontaktów atomowych i kryteria odległościowo-kątowe. Własny detektor C++ powinien być testowany różnicowo względem obu narzędzi, zamiast być oceniany wyłącznie „na oko”. ([OUP Academic][5])

Raport `bias_site.json` powinien zawierać co najmniej:

```json
{
  "source": "oracle_native",
  "feature_type": "ligand_acceptor",
  "target_residue": "ASP_123",
  "target_atom": "OD1",
  "ligand_atom_map": 17,
  "center": [12.34, 45.67, 8.91],
  "donor_acceptor_distance_A": 2.91,
  "angle_deg": 161.0,
  "occupancy": 1.0,
  "altloc": null,
  "selection_rank": 1,
  "selection_rule_version": "hbond-v1"
}
```

Nie należy wybierać miejsca na podstawie tego, które daje najlepszy docking. To byłoby strojenie na zbiorze testowym.

## Bias receptorowy bez badanego liganda

Docelowa wersja powinna generować pole z geometrii receptora:

* dla receptora-donora: miejsce dla ligandowego akceptora wzdłuż wektora D–H,
* dla receptora-akceptora: jedno lub kilka miejsc ligandowego donora zgodnych z geometrią wolnych par,
* dla pierścienia aromatycznego: centroid, wektor normalny i preferowana odległość,
* dla karboksylanu: dwa równoważne atomy jako grupa logiczna OR,
* dla metalu: geometria koordynacyjna i dozwolone typy ligandowych atomów,
* dla reszty o niepewnej protonacji: kilka alternatywnych hipotez z wagami.

„Dokowanie do konkretnego aminokwasu” powinno więc oznaczać dokowanie do zestawu konkretnych cech interakcyjnych tej reszty, a nie do abstrakcyjnego środka reszty.

---

# 8. Wody, jony i pozostałe HETATM

## Wariant podstawowy: receptor dry

Dla pierwszego benchmarku zastosowałbym jednoznaczny protokół:

* usunąć wszystkie wody krystalograficzne,
* usunąć bufor, glicerol, DMSO, siarczany i inne dodatki krystalizacyjne,
* usunąć ligand referencyjny,
* zachować kompletne wybrane białko,
* zachować wyłącznie funkcjonalne metale, kofaktory lub jony wymagane przez mechanizm kieszeni.

Nie stosuj reguły „zostawiamy wszystkie jony”. `SO4` jest jonem, ale często jest składnikiem warunków krystalizacji, a nie funkcjonalnym elementem receptora.

W 6JQR poza ligandem C6F występują m.in. `SO4`, `GOL` i `CXS`; w 5N9R występują `SO4`, `GOL` i `DMS`. Nie powinny automatycznie przechodzić do receptora dockingowego. ([ebi.ac.uk][6])

## Wariant dodatkowy: wet-strict

Nie ma jednej uniwersalnej reguły identyfikacji „potrzebnych” wód. Metody predykcji dobrze odnajdują część pozycji wód, ale znacznie gorzej przewidują ich wpływ termodynamiczny. ([Springer Nature Link][7])

Wody można dodać jako osobne badanie czułości, stosując wcześniej zamrożone kryteria, przykładowo:

* occupancy ≥ 0.8,
* brak altloc,
* B-factor niewyższy istotnie od lokalnego otoczenia białka,
* co najmniej dwa kontakty z białkiem,
* mała ekspozycja na rozpuszczalnik,
* konserwacja w co najmniej dwóch niezależnych strukturach tego receptora,
* lokalizacja zgodna po superpozycji struktur.

Woda wybrana dlatego, że mostkuje dokładnie badany ligand krystalograficzny, jest kolejnym wariantem oracle i tak trzeba ją oznaczyć.

AutoDock Vina ma osobny protokół hydrated docking wykorzystujący modyfikowane mapy AD4 i jawne, wypieralne wody związane z ligandem. Nie należy mieszać tego protokołu z podstawowym benchmarkiem dry bez osobnego ramienia eksperymentalnego. ([autodock-vina.readthedocs.io][8])

---

# 9. Jedna jednostka, łańcuchy i assembly

Nie wybieraj zawsze „jednego łańcucha” wyłącznie dla ujednolicenia plików. Wybieraj minimalną kompletną jednostkę biologiczną zawierającą całą kieszeń.

Reguła:

1. Sprawdzić preferowany biological assembly.
2. Sprawdzić, czy atomy innego łańcucha uczestniczą w kieszeni.
3. Jeżeli kieszeń jest wewnątrz jednego monomeru, wybrać jedną kopię.
4. Jeżeli kieszeń leży na interfejsie, zachować wszystkie wymagane łańcuchy.
5. Kopii z tej samej komórki krystalicznej nie traktować jako niezależnych kompleksów.

Dla Twojego zestawu:

* **6JQR**: jedna kopia FLT3, łańcuch A, preferowany assembly monomeryczny. ([ebi.ac.uk][6])
* **3CS9**: cztery kopie ABL1, łańcuchy A–D, ale preferowane assemblies są monomeryczne. Ligand NIL w łańcuchu A ma najlepszy ranking dopasowania do gęstości spośród czterech kopii, więc A jest rozsądnym prerejestrowanym wyborem. ([ebi.ac.uk][9])
* **5N9R**: dwie kopie USP7, A i B, preferowany assembly monomeryczny. Ligand 8RN w B ma minimalnie lepszy ogólny ranking dopasowania niż w A, więc przy regule „najlepszy ligand fit” wybrałbym B. ([ebi.ac.uk][10])

6JQR wymaga dodatkowej ostrożności. Ligand C6F ma RSR 0.20, RSCC 0.87 i jedynie 18. percentyl rankingu dopasowania. To nie dyskwalifikuje struktury, ale czyni ją raczej testem trudnym niż najlepszym przypadkiem do walidacji wzorcowej. Rozważyłbym oznaczenie jej jako `stress_target` albo zastąpienie w benchmarku potwierdzającym strukturą o lepszej jakości liganda. ([RCSB PDB][11])

---

# 10. Czym jest luka w PDB

Luka oznacza zwykle, że określonych reszt albo atomów nie udało się wymodelować z danych eksperymentalnych. Nie oznacza automatycznie luki w kieszeni.

mmCIF ma osobne kategorie dla reszt nieobserwowanych lub o zerowej occupancy. ([mmcif.k8s.wwpdb.org][12])

Proponuję automatyczny quality gate:

* brakujący atom ciężki reszty znajdującej się ≤ 6 Å od liganda referencyjnego — `FAIL` lub ręczny przegląd,
* brakujący fragment łańcucha, którego reszty flankujące znajdują się ≤ 8 Å od kieszeni — `FAIL`,
* luka daleko od kieszeni — raportować, ale nie wykluczać,
* altloc w kieszeni — wybrać deterministycznie najwyższą occupancy i zgłosić,
* zerowa occupancy w kieszeni — wykluczyć strukturę z głównego benchmarku.

Nie odbudowywałbym brakującej pętli lub side-chainu w głównym wariancie. Modelowanie dodaje nową warstwę niepewności. Lepsze są dwa osobne ramiona:

* `experimental_coordinates_only`,
* `repaired_structure_sensitivity`.

---

# 11. Przygotowanie receptora i liganda

## Receptor

Meeko powinno być główną współczesną ścieżką przygotowania. Legacy `prepare_receptor4.py` i `prepare_gpf4.py` warto zachować wyłącznie jako ścieżkę kompatybilności.

Nie porównuj dwóch metod przygotowania receptora w tym samym podstawowym eksperymencie bias versus no-bias. Wtedy nie byłoby wiadomo, czy wynik zmienił bias, czy preparacja.

Główna ścieżka:

```text
mmCIF
→ wybór assembly/łańcucha
→ kontrola missing atoms/altloc/occupancy
→ usunięcie dodatków
→ protonacja i przegląd sieci H-bond
→ Meeko receptor PDBQT
→ GPF
→ AutoGrid4 maps
```

Meeko obecnie odpowiada za typowanie atomów, ładunki, wiązania rotowalne i przygotowanie PDBQT dla Vina/AutoDock-GPU. ([meeko.readthedocs.io][13])

## Ligand

W Twoim zapisie jest ważny błąd:

> `RDKit AddHs (pH~7.4)`

`AddHs` jedynie dodaje wodory do istniejącego grafu. Nie przypisuje stanu protonacji przy pH 7.4 i nie wybiera tautomeru. ([rdkit.org][14])

Prawidłowy proces:

```text
SDF/CCD z bond orders i stereochemią
→ enumeracja protomerów/tautomerów
→ wybór lub zważenie mikrostanów dla zadanego pH
→ AddHs
→ niezależne wygenerowanie konformera 3D
→ Meeko
→ PDBQT
```

Oficjalny tutorial Vina wskazuje m.in. MolScrub do protonacji i enumeracji stanów, natomiast Meeko oczekuje SDF z jawnymi wodorami i współrzędnymi 3D. ([autodock-vina.readthedocs.io][15])

Rozdziel dwa pliki:

```text
ligand_reference.sdf
```

* niezmienne współrzędne krystalograficzne,
* wyłącznie do RMSD i opisu interakcji oracle;

```text
ligand_input.sdf
```

* mikrostan wybrany niezależnie,
* nowo wygenerowana konformacja,
* wejście do dockingu.

Nie należy przekazywać natywnej konformacji jako startowej, jeżeli którykolwiek etap może z niej korzystać. Vina wykonuje niezależne stochastic runs rozpoczynające się od losowych konformacji, ale nadal warto formalnie oddzielić referencję od wejścia. ([autodock-vina.readthedocs.io][2])

Do eksportu pozy z PDBQT użyj Meeko, ponieważ zachowuje informacje pozwalające odtworzyć bond orders i formalne ładunki z pierwotnego liganda. Oficjalna dokumentacja ostrzega, że zgadywanie ich na podstawie samego PDBQT może być niemożliwe. ([autodock-vina.readthedocs.io][15])

---

# 12. Box i „receptor ligand ± 15 Å”

Nie przycinaj receptora kolejno do 15, 10, 5 i 3 Å w podstawowym benchmarku.

Przycinanie receptora:

* zmienia środowisko elektrostatyczne,
* tworzy sztuczne granice,
* może usuwać reszty stabilizujące lokalną geometrię,
* dodaje kolejny czynnik eksperymentalny.

W głównym wariancie zachowaj cały wybrany domenowy receptor. Ograniczaj przestrzeń przeszukiwania przez box.

Prostopadłościan jest właściwszy niż sześcian, ponieważ nie marnuje ewaluacji w pustej przestrzeni.

Dla oracle redockingu:

1. Wyznaczyć bounding box ciężkich atomów natywnego liganda.
2. Dodać stały padding per oś, np. 5 Å z każdej strony.
3. Dopasować rozmiary do spacingu AutoGrid.
4. Sprawdzić, czy ligand i wszystkie miejsca biasu mają bezpieczny margines od granicy.
5. Zastosować dokładnie ten sam box we wszystkich warunkach.

Dla testu generalizacji box nie może być wyprowadzany z badanego liganda. Powinien pochodzić z:

* innych ligandów,
* listy reszt kieszeni,
* niezależnego pocket detectora,
* albo znanego miejsca funkcjonalnego receptora.

Oficjalna dokumentacja Vina potwierdza, że dla AD4 wymagane są mapy AutoGrid i wywołanie `--scoring ad4`; podaje również spacing 0.375 Å w przykładzie przygotowania map. ([autodock-vina.readthedocs.io][15])

---

# 13. Parametry biasu

## Wersja zgodnościowa

Dla pilota możesz zamrozić:

```yaml
potential: adbias_legacy_gaussian
Vset_kcal_mol: -1.5
radius_A: 1.0
```

ale opisz to jako:

> parametr kompatybilności z istniejącym przykładem,

a nie jako wartość wynikającą z „energii najsłabszego wiązania wodorowego”.

Scoring AD4/Vina jest modelem przybliżonym, a bias jest prior-em sterującym przeszukiwaniem. `Vset` nie jest pomiarem energii konkretnego fizycznego wiązania wodorowego.

Instrukcja AutoDock Bias opisuje około `-2 kcal/mol` jako silny bias i podaje orientacyjny zakres promienia około 0.6–1.2 Å. ([AutoDock Bias][16])

## Strojenie docelowe

Parametry stroić wyłącznie na oddzielnym zbiorze development:

```text
Vset ∈ {-0.5, -1.0, -1.5, -2.0, -3.0} kcal/mol
r    ∈ {0.6, 0.8, 1.0, 1.2, 1.5} Å
```

Następnie zamrozić jeden globalny zestaw na typ interakcji:

* osobny dla H-bond,
* osobny dla aromatic,
* osobny dla hydrophobic,
* osobny dla metal coordination.

Nie stosować innych parametrów dla każdego receptora, chyba że jest to jawnie badany model adaptacyjny uczony wyłącznie na danych treningowych.

Przy wielu miejscach wprowadź maksymalny budżet nagrody. Inaczej receptor z sześcioma miejscami otrzyma silniejszą ingerencję w funkcję celu niż receptor z jednym miejscem.

---

# 14. Lepsza matematyczna reprezentacja biasu

Legacy AdBias jest izotropowy i zależy tylko od odległości. Dla wiązania wodorowego to uproszczenie, ponieważ interakcja jest kierunkowa.

W C++ wprowadziłbym ogólną funkcję:

[
b_j(x,u)=V_j
\exp\left[
-\frac12(x-c_j)^T\Sigma_j^{-1}(x-c_j)
\right]
a_j(u)
]

gdzie:

* (x) — pozycja cechy ligandowej,
* (c_j) — środek miejsca,
* (\Sigma_j) — macierz opisująca niepewność i anizotropię,
* (u) — orientacja cechy,
* (a_j(u)\in[0,1]) — zgodność kątowa,
* (V_j<0) — maksymalna nagroda.

Dzięki temu:

* H-bond może mieć wąski rozkład poprzeczny i szerszy wzdłuż osi,
* aromatic stacking może uwzględniać centroid i normalną pierścienia,
* metal coordination może uwzględniać geometrię koordynacyjną,
* niepewne miejsce może mieć większą kowariancję zamiast arbitralnie dużego promienia.

## Alternatywne miejsca: logika OR

Dwa równoważne atomy karboksylanu albo kilka możliwych miejsc nie powinny sumować nagród. To sztucznie premiowałoby liczbę hipotez.

Dla alternatyw użyj gładkiego minimum:

[
B_{\mathrm{OR}}=
-\tau\log\left[
\frac1M\sum_{j=1}^{M}
\exp\left(-\frac{b_j}{\tau}\right)
\right]
]

## Interakcje wymagane jednocześnie: AND

Dla interakcji obowiązkowych można sumować wkłady, ale z limitem:

[
B_{\mathrm{AND}}=
\max\left(B_{\min},\sum_j b_j\right)
]

## Dalsze funkcjonalności

Warto obsłużyć:

* `OR`,
* `AND`,
* `k-of-n`,
* maksymalny budżet energii,
* dodatnie exclusion volumes,
* harmonogram annealingu biasu,
* confidence weights,
* grupy alternatywnych protonacji,
* propagację niepewności miejsca,
* source provenance.

rDock już reprezentuje restraints przez pozycję, typ cechy i promień tolerancji. Nowość Waszego rozwiązania powinna iść dalej: kierunkowość, niepewność, logika grup, dekompozycja energii i integracja z przeszukiwaniem Vina. ([rdock.github.io][17])

---

# 15. Porównanie AD4 i Vina

Trzeba rozdzielić trzy rzeczy:

1. algorytm przeszukiwania,
2. bazową funkcję scoringową,
3. sposób implementacji biasu.

`vina --scoring ad4` oznacza:

* algorytm przeszukiwania Vina,
* funkcję AD4 opartą na mapach AutoGrid.

Nie oznacza użycia klasycznego algorytmu AutoDock4 LGA. Oficjalna dokumentacja Vina przedstawia AD4 i Vina jako alternatywne force fields w tym samym silniku Vina. ([autodock-vina.readthedocs.io][15])

Proponowana macierz:

| Wariant | Search        | Base scoring | Bias                  |
| ------- | ------------- | ------------ | --------------------- |
| A       | Vina          | AD4 maps     | modyfikacja map       |
| B       | Vina          | Vina         | bezpośredni człon C++ |
| C       | AutoDock4 LGA | AD4 maps     | modyfikacja map       |
| D       | AutoDock-GPU  | AD4 maps     | modyfikacja map       |

Kolejność prac:

1. Wariant A — zgodność ze starym AdBias.
2. Wariant B — właściwa nowa funkcjonalność.
3. Wariant C — oddzielenie efektu search algorithm.
4. Wariant D — wydajność i duża skala.

Nie porównuj silników przy tym samym `exhaustiveness`, ponieważ parametr nie ma wspólnej interpretacji między różnymi algorytmami. Porównuj:

* accuracy przy kilku budżetach czasu,
* accuracy względem liczby ewaluacji,
* krzywe success–CPU time.

Na dzień 22 lipca 2026 oficjalna strona wydań oznacza Vina 1.2.7 jako najnowsze wydanie; zawiera ono m.in. poprawkę `set_weights` dla AD4. Wcześniejsze 1.2.3 naprawiało również błąd wymiarów siatek AD4. Dlatego wersja i commit muszą być zapisane w manifeście. ([GitHub][18])

---

# 16. Minimalny pilot, który ma sens

Dla obecnych trzech PDB proponuję:

## Struktury

```yaml
targets:
  - pdb_id: 6JQR
    chain: A
    ligand: C6F
    role: stress_target

  - pdb_id: 3CS9
    chain: A
    ligand: NIL
    role: primary

  - pdb_id: 5N9R
    chain: B
    ligand: 8RN
    role: primary
```

## Przygotowanie

* canonical input: mmCIF,
* pełna wybrana domena/łańcuch,
* receptor dry,
* usunięte wszystkie dodatki krystalizacyjne,
* zachowane tylko funkcjonalne metale/kofaktory,
* Meeko jako główna preparacja,
* legacy ADT tylko w testach zgodności,
* jeden zamrożony mikrostan liganda,
* osobne `reference.sdf` i `input.sdf`,
* box = ligand bounding box + 5 Å na każdą stronę,
* identyczne mapy bazowe we wszystkich warunkach.

## Bias

Dla każdego targetu:

* jedna bezpośrednia, dobrze określona interakcja H-bond,
* tylko jeżeli przechodzi quality gate,
* wybór automatyczny i zapisany przed dockingiem,
* `Vset=-1.5 kcal/mol`,
* `r=1.0 Å`,
* parametry identyczne dla wszystkich targetów.

Jeżeli któryś target nie ma jednoznacznej, jakościowej interakcji tego typu, należy wymienić target, a nie wymuszać wybór.

## Warunki

```yaml
conditions:
  - no_bias
  - correct_oracle_bias
  - matched_decoy_bias
```

`matched_decoy` powinien mieć:

* ten sam typ,
* tę samą wartość `Vset`,
* ten sam promień,
* podobną dostępność przestrzenną,
* podobny zakres bazowych wartości mapy,
* lokalizację w tej samej kieszeni,
* minimum 4 Å od prawidłowego miejsca,
* brak zgodności z natywną cechą liganda.

Losowy punkt w pustej przestrzeni jest zbyt łatwą kontrolą.

## Parametry uruchomienia

Jako prerejestrowany punkt startowy:

```yaml
vina:
  version: 1.2.7
  scoring: ad4
  exhaustiveness: 32
  num_modes: 20
  energy_range: 10
  cpu: 1
  seeds: 13
```

`cpu: 1` upraszcza analizę czasu i ogranicza wpływ planowania wątków. Dla benchmarku wydajności można później dodać osobny wariant wielowątkowy.

Oficjalna dokumentacja Vina opisuje `exhaustiveness` jako liczbę niezależnych stochastic runs i pokazuje, że zwiększenie go z domyślnych 8 do 32 może stabilizować trudniejsze przypadki. ([autodock-vina.readthedocs.io][15])

Łącznie:

[
3\times3\times13=117
]

przebiegów na silnik.

---

# 17. Docelowa architektura C++

Nie przepisywałbym od razu całej chemii przygotowania receptorów i ligandów. Najpierw należy zastąpić te elementy, w których C++ daje realną przewagę.

```text
src/
  io/
    mmcif_reader.*
    ad4_map.*
    pdbqt.*
    sdf_mapping.*

  structure/
    assembly_selector.*
    altloc_resolver.*
    gap_qc.*
    heterogen_policy.*
    pocket_definition.*

  chemistry/
    feature_typing.*
    interaction_detector.*
    atom_mapping.*
    microstate_manifest.*

  bias/
    bias_site.*
    legacy_gaussian.*
    anisotropic_potential.*
    angular_terms.*
    logic_groups.*
    decoy_generator.*
    perturbation_generator.*
    map_patcher.*

  engine/
    vina_ad4_adapter.*
    vina_native_adapter.*
    autodock4_adapter.*

  evaluation/
    symmetry_rmsd.*
    interaction_recovery.*
    pose_validity.*
    failure_classifier.*

  orchestration/
    run_manifest.*
    experiment_matrix.*
    cache.*
    provenance.*

  app/
    biasbench.cpp
```

Proponowane biblioteki:

* Gemmi do mmCIF i struktur,
* RDKit C++ do chemii liganda, mapowania i RMSD,
* Eigen do macierzy potencjałów anizotropowych,
* `yaml-cpp` oraz JSON Schema do konfiguracji,
* CLI11 do interfejsu,
* spdlog do logów,
* Catch2 lub GoogleTest do testów.

Model danych:

```cpp
struct BiasSite {
    FeatureType feature_type;
    Eigen::Vector3d center;
    Eigen::Matrix3d covariance;
    std::optional<Eigen::Vector3d> orientation;
    double depth_kcal_mol;
    double angular_sigma_rad;
    LogicGroupId logic_group;
    BiasSource source;
    double confidence;
    ResidueId target_residue;
    std::vector<AtomId> evidence_atoms;
};
```

`BiasSource` powinien jawnie rozróżniać:

```text
oracle_native
cross_ligand
receptor_geometry
homolog_structure
experimental
md_hotspot
water_prediction
manual
```

---

# 18. Funkcjonalności, które rzeczywiście będą lepsze niż AdBias

Najbardziej wartościowy zakres:

1. **Tryb pełnej kompatybilności BPF i map AD4.**
2. **Bezpośredni człon biasu w C++ Vina**, bez konieczności edycji plików map.
3. **Anizotropowe i kątowe oddziaływania.**
4. **Dekompozycja `base_score`, `bias_score`, `total_score`.**
5. **Bias tylko podczas search oraz czysty rescoring.**
6. **AND/OR/k-of-n dla alternatywnych interakcji.**
7. **Limit całkowitej nagrody.**
8. **Annealing**, np. silniejszy bias na początku i słabszy podczas końcowej lokalnej optymalizacji.
9. **Positive penalties i exclusion volumes.**
10. **Obsługa donorów, akceptorów, aromatyczności, salt bridges, metali, halogen bonds, cation–π i water bridges.**
11. **Confidence i niepewność przestrzenna.**
12. **Automatyczne generowanie dopasowanych kontroli negatywnych.**
13. **Źródło i provenance każdego miejsca.**
14. **Receptor ensembles i miejsca konserwowane między konformacjami.**
15. **Automatyczna klasyfikacja search failure versus ranking failure.**
16. **Wbudowany benchmark i generowanie danych tidy/Parquet.**

To jest dużo istotniejsza innowacja niż samo przepisanie Python 2 do C++.

---

# 19. Testy oprogramowania

## Golden tests

Porównanie z Python 2 AdBias punkt po punkcie:

```text
tests/golden_adbias/
  acc_single/
  don_single/
  aro_single/
  custom_map/
  overlapping_sites/
  boundary_site/
  off_grid_center/
  zero_depth/
```

## Property tests

Własności matematyczne:

* `Vset=0` daje mapę identyczną,
* wartość biasu dąży do zera wraz z odległością,
* dla izotropowego biasu wynik nie zależy od kierunku,
* kolejność miejsc w OR nie zmienia wyniku,
* cap energii nigdy nie jest przekraczany,
* brak NaN/Inf,
* header, origin, spacing i dimensions mapy pozostają niezmienione,
* miejsce poza mapą powoduje kontrolowany błąd.

## Gradient tests

Dla bezpośredniego biasu w Vina:

[
\frac{\partial B}{\partial x_i}
]

porównać z centralną różnicą skończoną dla losowych konfiguracji.

## Differential tests

* interaction detector C++ kontra PLIP,
* interaction detector C++ kontra Arpeggio,
* RMSD C++ kontra RDKit,
* map patcher C++ kontra legacy AdBias.

## RMSD fixtures

Obowiązkowe przykłady:

* pierścień fenylowy,
* karboksylan,
* nitro,
* dwa równoważne halogeny,
* symetryczny ligand,
* chiralność odwrócona,
* poprawna chemia, ale błędna pozycja przestrzenna.

## Integration tests

```text
mmCIF
→ receptor selection
→ receptor PDBQT
→ ligand PDBQT
→ GPF
→ AutoGrid4
→ bias
→ Vina
→ SDF
→ RMSD
→ results.parquet
```

## Testy statystyczne

Na danych syntetycznych:

* pod hipotezą zerową częstość fałszywych alarmów powinna być zgodna z nominalnym poziomem,
* przedziały ufności powinny mieć właściwe coverage,
* model powinien odzyskiwać znany zasymulowany efekt,
* błędne potraktowanie seedów jako niezależnych powinno być wykrywane przez test regresyjny analizy.

---

# 20. Reproducibility

Każdy run powinien mieć identyfikator:

[
\mathrm{run_id}=
\mathrm{hash}(
\mathrm{config}+
\mathrm{inputs}+
\mathrm{tool\ versions}+
\mathrm{git\ commit}
)
]

Manifest powinien zawierać:

* SHA-256 każdego wejścia,
* PDB ID i datę pobrania,
* wybrany assembly i chain,
* listę usuniętych heterogenów,
* listę zachowanych metali/kofaktorów,
* wszystkie decyzje protonacyjne,
* wersje Meeko, RDKit, AutoGrid i Vina,
* commit nowego kodu,
* compiler i flagi,
* system operacyjny i obraz kontenera,
* CPU,
* seed,
* box,
* wszystkie parametry biasu,
* ostrzeżenia QC.

Dane surowe powinny być niezmienne. Każdy etap tworzy nowy artefakt, a nie nadpisuje wcześniejszego.

Proponowany układ wyników:

```text
data/raw/<PDBID>/
data/derived/<PDBID>/
runs/<run_id>/<PDBID>/<engine>/<condition>/<seed>/
results/
  structures.parquet
  sites.parquet
  runs.parquet
  poses.parquet
  failures.parquet
report/
  figures/
  tables/
```

---

# 21. Wersja planu po proofreadzie

> **Cel projektu**
> Opracowanie i walidacja implementacji C++ miękkiego, receptorowo-specyficznego biasu interakcyjnego dla AutoDock Vina i funkcji AD4. Metoda ma zachowywać zgodność z mapowym modelem AutoDock Bias, ale rozszerzać go o bezpośrednią integrację z Vina, kierunkowe potencjały, niepewność przestrzenną, logiczne grupowanie alternatywnych interakcji, kontrolę całkowitej nagrody oraz rozdzielenie biasu używanego do przeszukiwania od końcowego scoringu.
>
> **Etap pilotażowy**
> Dla trzech nowych kompleksów krystalograficznych 6JQR, 3CS9 i 5N9R zostanie przeprowadzony kompletny pipeline od danych mmCIF i SDF do przygotowania receptora, przygotowania liganda, wygenerowania map AutoGrid4, automatycznego wyznaczenia miejsca biasu, dockingu i oceny pozy.
>
> **Porównywane warunki**
> Dla każdego kompleksu zostaną wykonane trzy warunki: docking bez biasu, docking z prawidłowym biasem oracle oraz docking z dopasowanym błędnym biasem. Wszystkie pozostałe elementy pipeline’u, w tym przygotowanie receptora i liganda, box, mapy bazowe, parametry Vina, budżet obliczeniowy i lista seedów, pozostaną identyczne.
>
> **Pierwszorzędowy punkt końcowy**
> Pierwszorzędowym wynikiem będzie prawdopodobieństwo uzyskania fizycznie poprawnej pozy top-1 o ciężkoatomowym, symetrycznie skorygowanym RMSD ≤ 2 Å względem pozy krystalograficznej.
>
> **Drugorzędowe punkty końcowe**
> Zostaną ocenione: ciągły RMSD, success@1/2/3 Å, best-of-K, pozycja pierwszej natywopodobnej konformacji w rankingu, odzyskanie wskazanej interakcji, fizyczna poprawność pozy, czysty score AD4/Vina, wkład biasu, różnorodność pozy oraz koszt obliczeniowy.
>
> **Interpretacja**
> Trzy kompleksy będą traktowane jako pilot techniczny i mechanistyczny. Seedy będą replikami technicznymi zagnieżdżonymi w kompleksach. Wnioski dotyczące generalizacji zostaną zweryfikowane w osobnym benchmarku leave-one-chemotype-out z biasami wyprowadzonymi bez użycia badanego liganda.

Najbliższym krokiem powinno być zamrożenie `pilot.yaml`, zasad wyboru łańcucha, protokołu usuwania HETATM, algorytmu wyboru biasu, generatora matched decoy oraz planu analizy przed uruchomieniem pierwszego porównawczego dockingu.

[1]: https://academic.oup.com/bioinformatics/article/35/19/3836/5368528 "https://academic.oup.com/bioinformatics/article/35/19/3836/5368528"
[2]: https://autodock-vina.readthedocs.io/en/latest/faq.html "https://autodock-vina.readthedocs.io/en/latest/faq.html"
[3]: https://pubmed.ncbi.nlm.nih.gov/38425520/ "https://pubmed.ncbi.nlm.nih.gov/38425520/"
[4]: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html "https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html"
[5]: https://academic.oup.com/nar/article/43/W1/W443/2467865 "https://academic.oup.com/nar/article/43/W1/W443/2467865"
[6]: https://www.ebi.ac.uk/pdbe/entry/pdb/6jqr/index "https://www.ebi.ac.uk/pdbe/entry/pdb/6jqr/index"
[7]: https://link.springer.com/article/10.1007/s10822-019-00187-y "https://link.springer.com/article/10.1007/s10822-019-00187-y"
[8]: https://autodock-vina.readthedocs.io/en/latest/docking_hydrated.html "https://autodock-vina.readthedocs.io/en/latest/docking_hydrated.html"
[9]: https://www.ebi.ac.uk/pdbe/entry/pdb/3cs9/index "https://www.ebi.ac.uk/pdbe/entry/pdb/3cs9/index"
[10]: https://www.ebi.ac.uk/pdbe/entry/pdb/5n9r/index "https://www.ebi.ac.uk/pdbe/entry/pdb/5n9r/index"
[11]: https://www.rcsb.org/ligand-validation/6JQR/C6F "https://www.rcsb.org/ligand-validation/6JQR/C6F"
[12]: https://mmcif.k8s.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/pdbx_unobs_or_zero_occ_residues.html "https://mmcif.k8s.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/pdbx_unobs_or_zero_occ_residues.html"
[13]: https://meeko.readthedocs.io/ "https://meeko.readthedocs.io/"
[14]: https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html "https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html"
[15]: https://autodock-vina.readthedocs.io/en/latest/docking_basic.html "https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"
[16]: https://autodockbias.wordpress.com/wp-content/uploads/2019/03/adbias_user_guide.pdf "https://autodockbias.wordpress.com/wp-content/uploads/2019/03/adbias_user_guide.pdf"
[17]: https://rdock.github.io/documentation/html_docs/user-guide/pharmacophoric-restraints.html "https://rdock.github.io/documentation/html_docs/user-guide/pharmacophoric-restraints.html"
[18]: https://github.com/ccsb-scripps/AutoDock-Vina/releases "https://github.com/ccsb-scripps/AutoDock-Vina/releases"
