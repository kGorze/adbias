1. Tryb pełnej kompatybilności BPF i map AD4.
2. Bezpośredni człon biasu w C++ Vina, bez konieczności edycji plików map.
3. Anizotropowe i kątowe oddziaływania.
4. Dekompozycja base_score, bias_score, total_score.
5. Bias tylko podczas search oraz czysty rescoring.
6. AND/OR/k-of-n dla alternatywnych interakcji.
7. Limit całkowitej nagrody.
8. Annealing, np. silniejszy bias na początku i słabszy podczas końcowej lokalnej optymalizacji.
9. Positive penalties i exclusion volumes.
10. Obsługa donorów, akceptorów, aromatyczności, salt bridges, metali, halogen bonds, cation–π i water bridges.
11. Confidence i niepewność przestrzenna.
12. Automatyczne generowanie dopasowanych kontroli negatywnych.
13. Źródło i provenance każdego miejsca.
14. Receptor ensembles i miejsca konserwowane między konformacjami.
15. Automatyczna klasyfikacja search failure versus ranking failure.
16. Wbudowany benchmark i testy.


testy
1: Test 1. Centrum zgodne z siatką
2: Test 2. Centrum niezgodne z siatką
3: Test 3. Symetria radialna
4: Test 4. Punkt d=r
5: Test 5. Próg
6: Test 6. Dyskretna siatka wyrównana
7: Test 7. Granica mapy
8: Test 8. Nakładające się biasy
9: Test 9. PFP
10: Test 10. Porównanie algorytmów
