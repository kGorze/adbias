bias.tcl — ładuje receptor oraz bias jako dwa osobne molecule.
bias_bias.tcl — samodzielny loader biasu do dowolnego otwartego systemu.
bias_scene.pdb — zawiera wyłącznie 1132 pomocnicze atomy biasu, bez receptora.


mol new inny_receptor.pdb
source bias_bias.tcl

python3 -m benchmark.versions.new.bias_vizualization one \
  --map benchmark/versions/new/results/3CS9/receptor.A.map \
  --bias-file benchmark/versions/new/results/3CS9/bias.bpf \
  --bias-number 1 \
  --receptor benchmark/versions/new/results/3CS9/receptor_prepared.pdb \
  --output benchmark/versions/new/results/3CS9/bias_001.tcl

  vmd -e benchmark/versions/new/results/3CS9/bias_001.tcl
  source /home/kgorzelanczyk/adbias/benchmark/versions/new/results/3CS9/bias_001_bias.tcl


  python3 -m benchmark.versions.new.bias_vizualization all \
  --results-dir benchmark/versions/new/results