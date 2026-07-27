## KROK 0 - Uruchomienie generatora

```
python3 -m bias_vizualization all --systems 6JQR
```

ma trzy pliki
```
results/6JQR/receptor.A.map
results/6JQR/bias.bpf
results/6JQR/receptor_prepared.pdb
```

tworzone są
```
results/6JQR/bias_visualizations/
- bias_001_acc_scene.pdb
- bias_001_acc_bias.tcl
- bias_001_acc.tcl
```

## PLIK 1 - receptor.A.map
Generator czyta tylko pierwsze sześć linii, a konkretnie:
```
SPACING 0.375
NELEMENTS 58 98 70
CENTER -27.848 -11.124 -27.881
```
Ich znaczenie:

```
SPACING   = odstęp pomiędzy punktami siatki: 0.375 Å
NELEMENTS = liczba przedziałów siatki: 58 × 98 × 70
CENTER    = środek całej siatki AutoDock
```

Na tej podstawie Python potrafi wyznaczyć współrzędne każdego punktu siatki.

## PLIK 2 - bias.bpf
Cały plik w tym przypadku ma tylko dwie linie:
```
x y z Vset r type
-27.632 -5.511 -32.244 -1.50 1.00 acc
```

Znaczenie tych danych:
```
x, y, z = środek biasu
Vset    = energia biasu w środku
r       = promień definiujący bias
type    = rodzaj biasu, tutaj acceptor
```

## PLIK 3 - receptor_prepared.pdb
Fragment:
```
ATOM      1  C   TYR A 572     -32.314   7.731 -17.275                       C
ATOM      2  O   TYR A 572     -32.467   8.959 -17.267                       O
ATOM      3  CA  TYR A 572     -33.499   6.791 -17.016                       C
ATOM      4  HA  TYR A 572     -33.121   5.829 -16.638                       H
ATOM      5  N   TYR A 572     -34.389   7.322 -15.965                       N
ATOM      6  H   TYR A 572     -34.002   7.091 -15.039                       H
ATOM      7  CB  TYR A 572     -34.250   6.526 -18.347                       C
ATOM      8 HB2  TYR A 572     -35.024   5.765 -18.168                       H
ATOM      9 HB3  TYR A 572     -34.684   7.481 -18.678                       H
ATOM     10  CG  TYR A 572     -33.329   6.032 -19.458                       C
ATOM     11 CD1  TYR A 572     -32.826   4.721 -19.444                       C
ATOM     12 HD1  TYR A 572     -33.127   4.037 -18.637                       H
...
```
Generator nie zmienia tego pliku. Zapisuje jego ścieżkę w głównym skrypcie Tcl, aby VMD mógł wczytać receptor.

## KROK 4 - python liczy scene
Na podstawie receptor.A.map i bias.bpf Python wyznacza między innymi:
```
środek biasu:
(-27.632, -5.511, -32.244)

najbliższy punkt siatki:
(-27.473, -5.499, -32.381)

promień powierzchni 1/e:
1.00

promień powierzchni epsilon:
około 2.24

kolejne zaakceptowane punkty siatki:
(-28.223, -6.249, -34.256)
(-27.848, -6.249, -34.256)
(-27.473, -6.249, -34.256)
...
```

W pamięci programu powstają obiekty zbliżone do:
```
Sphere(
    center=(-27.632, -5.511, -32.244),
    radius=0.1,
    color="purple",
    group="bias_center",
)

Sphere(
    center=(-27.473, -5.499, -32.381),
    radius=0.08,
    color="cyan",
    group="nearest_grid_point",
)

Sphere(
    center=(-27.632, -5.511, -32.244),
    radius=1.0,
    color="green",
    group="one_over_e_bias_surface",
)
```
To nadal nie są obiekty VMD. Jest to niezależny od VMD opis sceny w Pythonie.

## PLIK 4 - bias_001_acc_scene.pdb
Python zamienia obiekty sceny na pomocniczy PDB:
```
HETATM    1  V   BCT Z   1     -27.632  -5.511 -32.244  1.00  0.10      BVIZ V
HETATM    2  V   NGP Z   2     -27.473  -5.499 -32.381  1.00  0.08      BVIZ V
HETATM    3  V   BRS Z   3     -27.632  -5.511 -32.244  1.00  1.00      BVIZ V
HETATM    4  V   EPS Z   4     -27.632  -5.511 -32.244  1.00  2.24      BVIZ V
HETATM    5  V   GSP Z   5     -27.473  -5.499 -32.381  1.00  0.05      BVIZ V
HETATM    6  V   GSP Z   5     -27.098  -5.499 -32.381  1.00  0.05      BVIZ V
HETATM    7  V   ACY Z   6     -28.223  -6.249 -34.256  1.00  0.04      BVIZ V
HETATM    8  V   ACY Z   6     -27.848  -6.249 -34.256  1.00  0.04      BVIZ V
HETATM    9  V   ACY Z   6     -27.473  -6.249 -34.256  1.00  0.04      BVIZ V
HETATM   10  V   ACY Z   6     -27.098  -6.249 -34.256  1.00  0.04      BVIZ V
...
HETATM  895  V   ACR Z   8     -27.098  -5.499 -32.006  1.00  0.04      BVIZ V
HETATM  896  V   ACR Z   8     -27.848  -5.124 -32.006  1.00  0.04      BVIZ V
HETATM  897  V   ACR Z   8     -27.473  -5.124 -32.006  1.00  0.04      BVIZ V
HETATM  898  V   ACR Z   8     -27.473  -5.499 -31.631  1.00  0.04      BVIZ V
CONECT    5    6
END
```
Poszczególne pola są tutaj wykorzystywane w specjalny sposób:
```
                 GRUPA BCT  |    wspolrzedne XYZ     | beta | segment BVIZ | 
HETATM    1  V   BCT Z   1     -27.632  -5.511 -32.244  1.00  0.10          BVIZ V
```

Skróty grup:
```
BCT = bias center
NGP = nearest grid point
BRS = powierzchnia biasu 1/e
EPS = powierzchnia epsilon
GSP = odcinek pokazujący jeden krok siatki
ACY = zaakceptowane punkty, kolor żółty
ACO = zaakceptowane punkty, kolor pomarańczowy
ACR = zaakceptowane punkty, kolor czerwony
REJ = odrzucone punkty
```


## PLIK 5 - bias_001_acc_bias.tcl
Ten plik mówi VMD, jak wyświetlić dane z pomocniczego PDB:
```
# Reusable bias overlay. Source this file after loading any system.
source "/home/kgorzelanczyk/adbias/benchmark/versions/new/bias_vizualization/renderer.tcl"

set bias_file "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/6JQR/bias_visualizations/bias_001_acc_scene.pdb"

set bias_molecule [
    ::bias_vizualization::load_bias \
        $bias_file \
        "6JQR_bias_001_acc_bias" \
        0.35
]

# Representation: Bias center
::bias_vizualization::add_representation \
    $bias_molecule \
    "Bias center" \
    BCT \
    bias_center \
    "VDW 1.0 20" \
    purple

# Representation: Nearest grid point
::bias_vizualization::add_representation \
    $bias_molecule \
    "Nearest grid point" \
    NGP \
    nearest_grid_point \
    "VDW 1.0 16" \
    cyan

# Representation: Bias 1/e isosurface
::bias_vizualization::add_representation \
    $bias_molecule \
    "Bias 1/e isosurface" \
    BRS \
    one_over_e_bias_surface \
    "VDW 1.0 30" \
    green

# Representation: Epsilon energy isosurface
::bias_vizualization::add_representation \
    $bias_molecule \
    "Epsilon energy isosurface" \
    EPS \
    epsilon_energy_surface \
    "VDW 1.0 30" \
    red

# Representation: One grid-spacing step
::bias_vizualization::add_representation \
    $bias_molecule \
    "One grid-spacing step" \
    GSP \
    grid_spacing \
    "Bonds 0.075 12.0" \
    white

mol top $bias_molecule
```

Rzeczywisty wygenerowany plik zapisuje każde wywołanie w jednej linii. Powyżej rozbiłem je dla czytelności.

Przykład:
```
BCT bias_center "VDW 1.0 20" purple
```
oznacza:
```
wybierz rekordy PDB z resname BCT
nazwij selekcję bias_center
pokaż je jako kule VDW
ustaw rozdzielczość kuli na 20
ustaw kolor purple
```

## PLIK 6 - renderer.tcl 
To wspólna biblioteka wykonywana wewnątrz VMD. Nie zawiera współrzędnych konkretnego biasu.

```
namespace eval ::bias_vizualization {
    variable material_name "BiasVisualization"
}

proc ::bias_vizualization::load_receptor {
    receptor_file molecule_name
} {
    if {![file isfile $receptor_file]} {
        error "Receptor PDB does not exist: $receptor_file"
    }

    set molecule [
        mol new $receptor_file type pdb waitfor all
    ]

    mol rename $molecule $molecule_name
    return $molecule
}

proc ::bias_vizualization::load_bias {
    bias_file molecule_name opacity
} {
    variable material_name

    set molecule [
        mol new $bias_file type pdb waitfor all autobonds off
    ]

    mol rename $molecule $molecule_name
    mol delrep 0 $molecule

    set visualization [
        atomselect $molecule "segname BVIZ"
    ]

    $visualization set radius [
        $visualization get beta
    ]

    $visualization delete

    if {[lsearch -exact [material list] $material_name] < 0} {
        material add $material_name
    }

    material change opacity $material_name $opacity
    return $molecule
}
```

Najważniejsza operacja:
```
$visualization set radius [$visualization get beta]
```
W pliku sceny było:
```
HETATM 1 ... BCT ... 1.00 0.10 ... BVIZ
```
VMD początkowo interpretuje 0.10 jako pole beta. Renderer kopiuje je do pola radius, dlatego reprezentacja VDW rysuje kulę o promieniu 0.10.

Druga procedura tworzy reprezentacje:
```
proc ::bias_vizualization::add_representation {
    molecule label pdb_resname selection_name style color_name
} {
    variable material_name

    set atoms [
        atomselect $molecule "resname $pdb_resname"
    ]

    $atoms set resname $selection_name
    $atoms delete

    color Resname $selection_name $color_name

    mol representation {*}$style
    mol color ResName
    mol selection "resname $selection_name"
    mol material $material_name
    mol addrep $molecule

    set representation [
        expr {[molinfo $molecule get numreps] - 1}
    ]

    puts "VMD representation $representation: $label"
}
```

Dla BCT procedura wykonuje logicznie:
```
atomselect $bias_molecule "resname BCT"
color Resname bias_center purple
mol representation VDW 1.0 20
mol selection "resname bias_center"
mol addrep $bias_molecule
```

## PLIK 7 - bias_001_acc.tcl

To główny plik, który uruchamiamy w VMD:
```
# Receptor and bias are loaded as separate VMD molecules.

source "/home/kgorzelanczyk/adbias/benchmark/versions/new/bias_vizualization/renderer.tcl"

set receptor_file "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/6JQR/receptor_prepared.pdb"

set receptor_molecule [
    ::bias_vizualization::load_receptor \
        $receptor_file \
        "6JQR_bias_001_acc_receptor"
]

source "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/6JQR/bias_visualizations/bias_001_acc_bias.tcl"

mol top $receptor_molecule
```

można uruchomić przez -e albo przez VMD normalnie jako source
```
vmd -e benchmark/versions/new/results/6JQR/bias_visualizations/bias_001_acc.tcl
```
lub
```
source "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/6JQR/bias_visualizations/bias_001_acc.tcl"
```

## Co dokładnie dzieje się po source bias_001_acc.tcl
kolejność wykonywania rzeczy
```
1. VMD otwiera bias_001_acc.tcl
2. bias_001_acc.tcl ładuje renderer.tcl
3. renderer.tcl definiuje:
   |-- load_receptor
   |-- load_bias
   `-- add_representation
4. load_receptor ładuje receptor_prepared.pdb
5. bias_001_acc.tcl ładuje bias_001_acc_bias.tcl
6. bias_001_acc_bias.tcl ponownie ładuje renderer.tcl
7. load_bias ładuje bias_001_acc_scene.pdb
8. beta z PDB zostaje skopiowana do radius
9. add_representation tworzy reprezentację BCT
10. add_representation tworzy reprezentację NGP
11. add_representation tworzy reprezentację BRS
12. add_representation tworzy reprezentację EPS
13. add_representation tworzy reprezentację GSP
14. add_representation tworzy reprezentacje ACY/ACO/ACR
15. VMD pokazuje receptor i nakładkę biasu jako dwie osobne molekuły
```