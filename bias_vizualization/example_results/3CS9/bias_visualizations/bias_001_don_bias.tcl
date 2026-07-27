# Reusable bias overlay. Source this file after loading any system.
source "/home/kgorzelanczyk/adbias/benchmark/versions/new/bias_vizualization/renderer.tcl"
set bias_file "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/3CS9/bias_visualizations/bias_001_don_scene.pdb"
set bias_molecule [::bias_vizualization::load_bias $bias_file "3CS9_bias_001_don_bias" 0.35]

# Representation: Bias center
::bias_vizualization::add_representation $bias_molecule "Bias center" BCT bias_center "VDW 1.0 20" purple
# Representation: Nearest grid point
::bias_vizualization::add_representation $bias_molecule "Nearest grid point" NGP nearest_grid_point "VDW 1.0 16" cyan
# Representation: Bias 1/e isosurface
::bias_vizualization::add_representation $bias_molecule "Bias 1/e isosurface" BRS one_over_e_bias_surface "VDW 1.0 30" green
# Representation: Epsilon energy isosurface
::bias_vizualization::add_representation $bias_molecule "Epsilon energy isosurface" EPS epsilon_energy_surface "VDW 1.0 30" red
# Representation: One grid-spacing step
::bias_vizualization::add_representation $bias_molecule "One grid-spacing step" GSP grid_spacing "Bonds 0.075 12.0" white
# Representation: Accepted points: low fraction
::bias_vizualization::add_representation $bias_molecule "Accepted points: low fraction" ACY accepted_low_fraction "VDW 1.0 4" yellow
# Representation: Accepted points: medium fraction
::bias_vizualization::add_representation $bias_molecule "Accepted points: medium fraction" ACO accepted_medium_fraction "VDW 1.0 4" orange
# Representation: Accepted points: high fraction
::bias_vizualization::add_representation $bias_molecule "Accepted points: high fraction" ACR accepted_high_fraction "VDW 1.0 4" red

mol top $bias_molecule
