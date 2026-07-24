# VMD backend. It only executes concrete drawing commands from generated Tcl.
namespace eval ::bias_vizualization {
    variable material_name "BiasVizMaterial"
}

proc ::bias_vizualization::prepare {receptor_file scene_name opacity} {
    variable material_name

    if {![file isfile $receptor_file]} {
        error "Receptor file does not exist: $receptor_file"
    }
    if {$opacity <= 0.0 || $opacity > 1.0} {
        error "Graphics opacity must be in the range (0, 1]"
    }

    set molecule [mol new $receptor_file type pdb waitfor all]
    mol rename $molecule $scene_name

    if {[lsearch -exact [material list] $material_name] < 0} {
        material add $material_name
    }
    material change opacity $material_name $opacity
    graphics $molecule material $material_name
    return $molecule
}

proc ::bias_vizualization::sphere {molecule center radius color resolution} {
    graphics $molecule color $color
    graphics $molecule sphere $center radius $radius resolution $resolution
}

proc ::bias_vizualization::point {molecule position color} {
    graphics $molecule color $color
    graphics $molecule point $position
}

proc ::bias_vizualization::line {molecule start end color width style} {
    graphics $molecule color $color
    graphics $molecule line $start $end width $width style $style
}

proc ::bias_vizualization::text {molecule position label size} {
    graphics $molecule color white
    graphics $molecule text $position $label size $size thickness 2
}
