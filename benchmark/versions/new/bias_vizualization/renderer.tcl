# VMD backend for separate receptor and reusable bias molecules.
namespace eval ::bias_vizualization {
    variable material_name "BiasVisualization"
}

proc ::bias_vizualization::load_receptor {receptor_file molecule_name} {
    if {![file isfile $receptor_file]} {
        error "Receptor PDB does not exist: $receptor_file"
    }

    set molecule [mol new $receptor_file type pdb waitfor all]
    mol rename $molecule $molecule_name
    return $molecule
}

proc ::bias_vizualization::load_bias {bias_file molecule_name opacity} {
    variable material_name

    if {![file isfile $bias_file]} {
        error "Bias PDB does not exist: $bias_file"
    }
    if {$opacity <= 0.0 || $opacity > 1.0} {
        error "Representation opacity must be in the range (0, 1]"
    }

    set molecule [mol new $bias_file type pdb waitfor all autobonds off]
    mol rename $molecule $molecule_name
    mol delrep 0 $molecule

    set visualization [atomselect $molecule "segname BVIZ"]
    if {[$visualization num] == 0} {
        $visualization delete
        error "Bias PDB does not contain the BVIZ segment"
    }
    $visualization set radius [$visualization get beta]
    $visualization delete

    if {[lsearch -exact [material list] $material_name] < 0} {
        material add $material_name
    }
    material change opacity $material_name $opacity
    return $molecule
}

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
