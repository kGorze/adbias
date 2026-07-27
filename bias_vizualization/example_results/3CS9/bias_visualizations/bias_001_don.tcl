# Receptor and bias are loaded as separate VMD molecules.
source "/home/kgorzelanczyk/adbias/benchmark/versions/new/bias_vizualization/renderer.tcl"
set receptor_file "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/3CS9/receptor_prepared.pdb"
set receptor_molecule [::bias_vizualization::load_receptor $receptor_file "3CS9_bias_001_don_receptor"]
source "/home/kgorzelanczyk/adbias/benchmark/versions/new/results/3CS9/bias_visualizations/bias_001_don_bias.tcl"
mol top $receptor_molecule
