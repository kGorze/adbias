# VMD backend for separate receptor and reusable bias molecules.


#Obecny przepływ wygląda tak:
#   1. Python oblicza geometrię biasu i buduje abstrakcyjną scenę: kule, punkty oraz linie.
#   2. Geometria zostaje zapisana jako pomocniczy PDB, np. bias_001_acc_scene.pdb.
#   3. Python generuje plik *_bias.tcl, opisujący reprezentacje VMD: selekcje, style i kolory.
#   4. Główny plik *.tcl ładuje receptor, a następnie nakładkę biasu.
#   5. Oba skrypty wykonują source renderer.tcl, aby uzyskać trzy procedury obsługujące VMD.


#mając takie coś nie zaśmiecamy przestrzeni nazw VMD
namespace eval ::bias_vizualization {
    variable material_name "BiasVisualization"
}

#TODO: do dokończenia jako poprawka, żeby te namespace były inne i można było wczytywać kilka biasó z różną przezroczystością, bo teraz wczytanie kolejnego biasu zmienia przezroczystość wszystkich biasów
namespace eval::bias_visualization_system {
    variable receptor_molecule
    variable bias_molecule
}

proc ::bias_vizualization::load_receptor {receptor_file molecule_name} {
    if {![file isfile $receptor_file]} {
        error "Receptor PDB does not exist: $receptor_file"
    }

    set molecule [mol new $receptor_file type pdb waitfor all]
    mol rename $molecule $molecule_name
    return $molecule
}

# wszystkie nakładki korzystają obecnie z jednego globalnego materiału BiasVisualization. załadowanie kolejnej nakładki z inną przezroczystością spowoduje zmianę przezroczystości wszystkich nakładek
proc ::bias_vizualization::load_bias {bias_file molecule_name opacity} {
    variable material_name

    if {![file isfile $bias_file]} {
        error "Bias PDB does not exist: $bias_file"
    }

    #TODO: sprawdzenie opacity, chyba powinno być w innym miejscu 
    if {$opacity <= 0.0 || $opacity > 1.0} {
        error "Representation opacity must be in the range (0, 1]"
    }

    set molecule [mol new $bias_file type pdb waitfor all autobonds off]
    mol rename $molecule $molecule_name
    mol delrep 0 $molecule

    #segment BVIZ identyfikacja dane wizualne
    set visualization [atomselect $molecule "segname BVIZ"]
    if {[$visualization num] == 0} {
        $visualization delete
        error "Bias PDB does not contain the BVIZ segment"
    }

    #beta jest jako promień kuli 
    $visualization set radius [$visualization get beta]
    $visualization delete

    if {[lsearch -exact [material list] $material_name] < 0} {
        material add $material_name
    }
    material change opacity $material_name $opacity
    return $molecule
}


# tworzenie selekcji i reprezetnacji dla danych z pythona
proc ::bias_vizualization::add_representation {
    molecule label pdb_resname selection_name style color_name
} {
    variable material_name

    set atoms [atomselect $molecule "resname $pdb_resname"]
    if {[$atoms num] == 0} {
        $atoms delete
        error "No bias atoms found for representation: $label"
    }
    $atoms set resname $selection_name
    $atoms delete

    color Resname $selection_name $color_name
    mol representation {*}$style
    mol color ResName
    mol selection "resname $selection_name"
    mol material $material_name
    mol addrep $molecule
    set representation [expr {[molinfo $molecule get numreps] - 1}]
    puts "VMD representation $representation: $label"
}
