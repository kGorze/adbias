from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem


analysis_dir = "/home/kgorzelanczyk/vs_dane/docking/analysis"
analysis_worst = "/home/kgorzelanczyk/vs_dane/docking/analysis/analysis_worst"
analysis_best = "/home/kgorzelanczyk/vs_dane/docking/analysis/analysis_best"

analysis_id = 0
#wczytanie danych z csv do programu
#format danych jakie są w csv
"""
ligand_id,smiles,score
0,CC(=O)N1CCC(O)C1,-4.318
1,Cc1ccnc(N)n1,-4.063
2,CN1CCNc2ccccc2C1,-4.595
"""
def import_molecules_from_csv(file_path):
    #data frame z csv
    df = pd.read_csv(file_path)
    molecules = []
    for index, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            molecules.append(mol)
    return molecules


# funkcja: zrobienie histogramu z wynikami dockingowymi i zapisanie do png
def plot_score_histogram(file_path, output_png, bins=20):
    df = pd.read_csv(file_path)
    plt.figure()
    plt.hist(df['score'], bins=bins)
    plt.xlabel("Wynik dokowania (kcal/mol)")
    plt.ylabel("Liczba ligandów")
    plt.title("Rozkład wyników dokowania")
    plt.savefig(output_png)
    plt.close()


# funkcja: filtorwanie po scoringu 20 systemow na plik results
def filter_top_scoring(file_path, output_path, top_n=20):
    df = pd.read_csv(file_path)
    top_df = df.sort_values("score").head(top_n)
    top_df.to_csv(output_path, index=False)
    return top_df

def filter_worst_scoring(file_path, output_path, worst_n=20):
    df = pd.read_csv(file_path)
    worst_df = df.sort_values("score", ascending=False).head(worst_n)
    worst_df.to_csv(output_path, index=False)
    return worst_df


# funkcja: do tych 20 systemow z najlepszym results zrobienie automatycznie skryptu tcl z wczytaniem receptora i ligandow zadokowanych
def generate_tcl_script(top_results, receptor_pdb, mol2_dir, output_tcl):
    lines = [f"mol new {receptor_pdb} type pdb"]
    for ligand_id in top_results['ligand_id']:
        mol2_path = f"{mol2_dir}/ligand_{ligand_id}_docked.mol2"
        lines.append(f"mol new {mol2_path} type mol2")

    Path(output_tcl).write_text("\n".join(lines) + "\n")


# definicja procedury liczącej kontakty i rysującej dla
# każdej pozy (klatki) osobny obiekt graficzny w VMD, plus synchronizacja z
# animacją. 
# 
# Kontakty są liczone przez  `measure contacts`, bo ligand i receptor są
# wczytane jako dwie różne molekuły
_TCL_CONTACTS_BODY = r'''

#wczytanie pakietu topotools
package require topotools

# ustawia aktywną klatkę molekuły liganda na klatkę 0, przed dalszymi operacjami.
molinfo $lig set frame 0

# atomselect tworzy selekcję atomów. frame now oznacza, że selekcja odnosi się do aktualnie ustawionej klatki
set _recsel [atomselect $rec all]
set _ligsel [atomselect $lig all frame now]

# funkcja selections2mol tworzy nową, osobną molekułę składającą się z atomów podanych selekcji.
set merged [::TopoTools::selections2mol [list $_recsel $_ligsel]]
mol rename $merged "receptor_ligand_merged"

#trzeba usuwać zasoby w VMD
$_recsel delete
$_ligsel delete

# selections2mol robi tylko SNAPSHOT jednej pozy, wiec scalona molekula bylaby
# statyczna.

# są dodane nowe klatki receptor stoi w miejscu, a atomy liganda
# dostaja wspolrzedne z kolejnych pozy. Atomy liganda leza w scalonej molekule
# za atomami receptora (bo laczylismy w kolejnosci [list $_recsel $_ligsel]).
set _nrec    [molinfo $rec get numatoms]
set _nlig    [molinfo $lig get numatoms]
set _nframes [molinfo $lig get numframes]
set _mlig [atomselect $merged "index $_nrec to [expr {$_nrec + $_nlig - 1}]"]

for {set f 1} {$f < $_nframes} {incr f} {
    # duplikuje klatke 0 (pelny receptor+ligand) i dokleja ja na koncu
    animate dup frame 0 $merged
    set _newframe [expr {[molinfo $merged get numframes] - 1}]

    # nadpisujemy wspolrzedne czesci ligandowej pozycja z klatki $f
    set _lsf [atomselect $lig all frame $f]
    $_mlig frame $_newframe
    $_mlig set {x y z} [$_lsf get {x y z}]
    $_lsf delete
}
$_mlig delete
puts "Scalona molekula ma [molinfo $merged get numframes] klatek (po jednej na poze)"

# Pobiera nazwę reszty (resname) liganda w sposób dynamiczny, zamiast wpisywać ją na sztywno. $_ls get resname zwraca listę resname dla wszystkich atomów (zwykle same powtórzenia tej samej wartości).
set _ls [atomselect $lig all]
set lig_resname [lindex [lsort -unique [$_ls get resname]] 0]
$_ls delete

# druga reprezentacja w scalonej molekule: kieszen wiazaca = wszystko < 5 A od liganda
mol addrep $merged
set repid [expr {[molinfo $merged get numreps] - 1}]


# selekcję (całe reszty mające choć jeden atom w promieniu 5 Å od liganda, czyli kieszeń wiążąca - "same residue as" zapobiega rysowaniu obciętych, wiszących w powietrzu fragmentów reszt), styl rysowania (Licorice z parametrami grubości i rozdzielczości) oraz schemat kolorowania (po nazwie atomu)
mol modselect $repid $merged "same residue as (all within 5 of resname $lig_resname)"
mol modstyle  $repid $merged Licorice 0.2 12 12
mol modcolor  $repid $merged Name

# KLUCZOWE: bez tego selekcja "within 5" jest liczona RAZ (dla klatki 0) i kieszen
# nie zmienia sie przy przewijaniu pozy. selupdate przelicza ja dla kazdej klatki.
mol selupdate $repid $merged on

# proc: dla jednej klatki ($frame) liganda znajduje reszty receptora w kontakcie
# i rysuje po jednej najkrótszej linii do każdej kontaktującej reszty.
proc make_pose_contacts {rec lig frame cutoff} {

    # Przełącza molekułę liganda na daną klatkę (pozę).
    molinfo $lig set frame $frame


    set ligsel [atomselect $lig "not hydrogen" frame $frame]
    set recsel [atomselect $rec "protein and not hydrogen"]

    set contacts [measure contacts $cutoff $recsel $ligsel]
    set rec_indices [lindex $contacts 0]
    set lig_indices [lindex $contacts 1]

    # pierwszy etap, który dla każdej reszty zapamiętujemy tylko najbliższą parę atomów
    set best_contacts [dict create]

    # foreach z dwiema zmiennymi i dwiema listami jednocześnie iteruje równolegle po obu listach. ri to kolejny indeks atomu receptora, li odpowiadający mu indeks atomu liganda z tej samej pary kontaktu.
    foreach ri $rec_indices li $lig_indices {

        set ra [atomselect $rec "index $ri"]
        set la [atomselect $lig "index $li" frame $frame]

        set receptor_xyz [lindex [$ra get {x y z}] 0]
        set ligand_xyz   [lindex [$la get {x y z}] 0]

        
        set receptor_info [lindex [$ra get {resname resid chain segname}] 0]

        #Dla pojedynczego atomu tworzy selekcję po indeksie, pobiera współrzędne XYZ oraz metadane reszty (nazwa, numer, łańcuch, segment). lassign to komenda Tcl rozdzielająca listę na osobne zmienne w jednej linii, odpowiednik rozpakowania krotki w Pythonie.
        lassign $receptor_info resname resid chain segname

        # vecsub i veclength to komendy wektorowe VMD. vecsub odejmuje wektory (współrzędne), veclength liczy długość wektora, czyli tutaj odległość euklidesową między atomem receptora a atomem liganda.
        set distance [veclength [vecsub $receptor_xyz $ligand_xyz]]


        #Buduje słownikdict exists sprawdza, czy klucz już istnieje w słowniku. 
        # dict get pobiera wartość dla klucza(odleglosc) Jeśli reszta jeszcze nie ma zapisanego kontaktu, albo nowy kontakt jest bliższy niż poprzedni, dict set nadpisuje wpis nową, bliższą parą. Efekt: na koniec pętli best_contacts zawiera dokładnie jeden, najbliższy kontakt na każdą kontaktującą resztę receptora. Na końcu każda tymczasowa selekcja atomu jest usuwana (delete), zgodnie z tą samą zasadą co wcześniej.
        set residue_key "$segname|$chain|$resid|$resname"

        if {![dict exists $best_contacts $residue_key]
            || $distance < [lindex [dict get $best_contacts $residue_key] 0]} {
            dict set best_contacts $residue_key \
                [list $distance $ligand_xyz $receptor_xyz $resname $resid $chain $segname]
        }

        $ra delete
        $la delete
    }

    $ligsel delete
    $recsel delete

    # mol new tworzy zupełnie nową, pustą molekułę służącą wyłącznie jako kontener na obiekty graficzne (linie, kule, tekst).
    set gfx [mol new]

    # mol rename nadaje jej czytelną nazwę z numerem pozy format "%02d" to formatowanie liczby na dwie cyfry z zerem wiodącym (np. 01, 02).
    mol rename $gfx [format "CONTACTS_pose_%02d" [expr {$frame + 1}]]
    
    # animate dup $gfx dodaje jedną klatkę do tej nowej, pustej molekuły (bez tego molekuła nie miałaby żadnej klatki i nie dałoby się na niej rysować).
    animate dup $gfx
    
    #zmiana koloru na żółty
    graphics $gfx color yellow

    # dane kontaktu
    dict for {residue_key data} $best_contacts {

        # Dla każdej najbliższej pary kontaktowej rysuje: cylinder (linię łączącą ligand z receptorem), dwie kule (na atomie receptora i liganda, różne promienie dla odróżnienia), oraz etykietę tekstową w połowie odległości (vecscale 0.5 [vecadd ...] liczy środek odcinka, suma wektorów razy 0.5) z nazwą reszty i odległością zaokrągloną do dwóch miejsc po przecinku.
        lassign $data distance ligand_xyz receptor_xyz resname resid chain segname

        # linia miedzy ligandem a receptorem
        graphics $gfx cylinder $ligand_xyz $receptor_xyz radius 0.10 resolution 12 filled yes
        
        # robienie kulek na atomach ligand i receptor
        graphics $gfx sphere $receptor_xyz radius 0.35 resolution 15
        graphics $gfx sphere $ligand_xyz radius 0.25 resolution 15

        # opis reszty i odległość
        set midpoint [vecscale 0.5 [vecadd $ligand_xyz $receptor_xyz]]

        if {$chain eq ""} {
            set label "${resname}${resid}"
        } else {
            set label "${resname}${resid}:${chain}"
        }

        # wypis tego
        graphics $gfx text $midpoint "$label [format %.2f $distance] A" size 0.8 thickness 1.0
    }

    puts "Poza [expr {$frame + 1}]: [dict size $best_contacts] kontaktujacych reszt"
    return $gfx
}

# array unset poseGraphics czyści całą tablicę poseGraphics, usuwając wszystkie jej elementy, gdyby istniały z poprzedniego uruchomienia skryptu. To zabezpieczenie przed powtórnym wykonaniem skryptu w tej samej sesji VMD.
array unset poseGraphics
set number_of_frames [molinfo $lig get numframes]
for {set frame 0} {$frame < $number_of_frames} {incr frame} {
    set poseGraphics($frame) [make_pose_contacts $rec $lig $frame $cutoff]

    # mol off ukrywa daną molekułę w widoku (wszystkie kontakty są od razu tworzone, ale na starcie niewidoczne).
    mol off $poseGraphics($frame)
}

# proc: pokazuje kontakty tylko dla aktualnie wyświetlanej klatki liganda.
proc update_pose_contacts {name1 name2 operation} {
    
    # global lig poseGraphics deklaruje, że procedura ma używać zmiennych globalnych lig i poseGraphics, a nie tworzyć lokalne kopie
    global lig poseGraphics merged
    set current_frame [molinfo $lig get frame]

    # scalona molekula ma tyle samo klatek co ligand — ustawiamy ja na te sama poze,
    # zeby merge i kieszen wiazaca szly w parze z animacja liganda
    if {[info exists merged] && $current_frame < [molinfo $merged get numframes]} {
        molinfo $merged set frame $current_frame
    }

    # Pętla wyłącza (mol off) wszystkie molekuły z kontaktami. Następnie info exists poseGraphics($current_frame) sprawdza, czy dla aktualnej klatki istnieje odpowiadający jej zestaw kontaktów
    foreach frame [array names poseGraphics] {
        mol off $poseGraphics($frame)
    }
    if {[info exists poseGraphics($current_frame)]} {
        mol on $poseGraphics($current_frame)
    }
}

# automatyczna synchronizację z suwakiem/animacją klatek liganda.
catch {
    # vmd_frame($lig) to specjalna zmienna globalna utrzymywana przez samo VMD, która automatycznie zmienia wartość, gdy użytkownik przesuwa suwak animacji
    trace remove variable vmd_frame($lig) write update_pose_contacts
}
trace add variable vmd_frame($lig) write update_pose_contacts

mol top $lig
molinfo $lig set frame 0

# Ustawia klatkę z powrotem na 0. Na końcu ręcznie wywołuje update_pose_contacts z trzema pustymi argumentami (spełniając wymaganą sygnaturę, mimo że nie są używane), żeby od razu na starcie pokazać kontakty dla pierwszej pozy, bez konieczności ruszania suwakiem.
update_pose_contacts {} {} {}
'''


# funkcja: zrobienie skryptu tcl, która tworzy w vmd wizualizacje z wczytaniem receptra i ligandów zadokowanych.
# każda pozycja ligandu ma być połączona z konkretnymi punktami. dla ligandu ma być połączenie z każdym aminokwasem, który może być w kontakcie z ligandem.
# kontakt jest zdefiniowany jako odległość mniejsza niż X do zdefiniowania w funkcji przez parametr Å między atomami ligandu a najbliższym atomem aminokwasu.
# dla kazdej klatki/pozycji ligandu ma być nowy obiekt w vmd, zeby mozna bylo sobie oglądać wszystkie w graphics jako nowa klatka a jednoczesnie dzilac z animacjami.
def generate_tcl_with_contacts(ligand_id, receptor_pdb, mol2_dir, output_tcl, contact_distance=4.0):
    mol2_path = f"{mol2_dir}/ligand_{ligand_id}_docked.mol2"

    # najpierw wczytujemy receptor potem ligand
    header = (
        f'set lig [mol new "{mol2_path}" type mol2 waitfor all]\n'
        f'mol rename $lig "ligand_{ligand_id}"\n'
        f"\n"
        f'set rec [mol new "{receptor_pdb}" type pdb waitfor all]\n'
        f'mol rename $rec "receptor"\n'
        f"\n"
        f"set cutoff {contact_distance}\n"
    )

    Path(output_tcl).write_text(header + _TCL_CONTACTS_BODY)
    return output_tcl


# funkcja: dla kazdego z 20 najlepszych ligandow tworzy podkatalog ligand_<id>/
# i zapisuje w nim powtarzalny skrypt vina_contacts_ligand_<id>.tcl.
def generate_tcl_for_20(results, receptor_pdb, mol2_dir, output_dir, contact_distance=4.0):
    generated = []
    for ligand_id in results['ligand_id']:
        ligand_dir = Path(output_dir) / f"ligand_{ligand_id}"
        ligand_dir.mkdir(parents=True, exist_ok=True)
        output_tcl = ligand_dir / f"vina_contacts_ligand_{ligand_id}.tcl"
        generate_tcl_with_contacts(
            ligand_id, receptor_pdb, mol2_dir, output_tcl, contact_distance=contact_distance
        )
        generated.append(output_tcl)
    return generated


if __name__ == "__main__":
    results_dir = "/home/kgorzelanczyk/vs_dane/docking/results"
    mol2_dir = f"{results_dir}/mol2"
    receptor_pdb = "/home/kgorzelanczyk/vs_dane/docking/receptor/6JQR_receptorFH.pdb"

    out_dir = f"{analysis_dir}/analysis_{analysis_id}"
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    results_csv = f"{results_dir}/results.csv"
    histogram_png = f"{out_dir}/scores_histogram.png"
    top20_csv = f"{out_dir}/top20_results.csv"
    worst20_csv = f"{out_dir}/worst20_results.csv"
    tcl_script = f"{out_dir}/load_top20.tcl"

    plot_score_histogram(results_csv, histogram_png)
    top20 = filter_top_scoring(results_csv, top20_csv, top_n=20)
    worst20 = filter_worst_scoring(results_csv, worst20_csv, worst_n=20)
    generate_tcl_script(top20, receptor_pdb, mol2_dir, tcl_script)
    generate_tcl_for_20(top20, receptor_pdb, mol2_dir, out_dir, contact_distance=4.0)
    generate_tcl_for_20(worst20, receptor_pdb, mol2_dir, out_dir, contact_distance=5.0)
