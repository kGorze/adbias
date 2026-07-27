# zrozumienie plików biasu
bias.tcl — ładuje receptor oraz bias jako dwa osobne molecule.
bias_bias.tcl — samodzielny loader biasu do dowolnego otwartego systemu.
bias_scene.pdb — zawiera wyłącznie 1132 pomocnicze atomy biasu, bez receptora.


# wczytywanie plików biasu
mol new inny_receptor.pdb
source bias_bias.tcl

# przykład użycia wizualizacji na systemach wykrywanych w moim przypadku 
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



  # interaction sites

  python interactions.py -i 3CS9 -c A -r 382



# wykrywanie systemów 
żeby to działało z wykrywaniem systemów, trzeba mieć je w results

/bias_visuzalization
/results
  - /system1
  - /system2
  - /system ...
  - /system n

# stałe które są wpisane w kodzie

1. 
generate_bias_visualization( ...
epsilon 0.01

)

2. 
generate_for_system(... 
epsilon: 0.01

)

3. 
generate_for_systems(... 
epsilon: 0.01

)

4. 
_corrected_axis_range():
tolerance = 1.0e-12

5. 
calculate_bias_geometry( ...
epsilon 0.01)

# discover systems

w disocver systems jest ustalone tak, że nazwy muszą mieć konkretny format:
```python
    map_filename: str = "receptor.A.map",
    bias_filename: str = "bias.bpf",
    receptor_filename: str = "receptor_prepared.pdb",
```

# czy można obecnie używać wielu biasów jednocześnie? 
nie, przez to, że wczytywanie biasów jest na jednej reprezentacji

# czy można robić wiele biasów na systemie?
nie sprawdzone, autorzy wskazywali, że tak jednak wtedy jest tylko jedno minimum a biasy nie są kumulatywne. nie działa wyłącznie przypadek używania biasu aromatycznego

# jak używać API?
używać tylko metod które są bez "_*", czyli wszystko co jest w API

```python
from .api import (
    SystemFiles,
    discover_systems,
    generate_bias_visualization,
    generate_for_system,
    generate_for_systems,
    parse_bias_file,
    renderer_path,
)
```

# instalacja

```bash
python -m pip install ./bias_vizualization
```


po zmianie plikow zrodlowych znowu trzeba zainstalowac paczke

```bash
python -m pip install --editable \
  ./benchmark/versions/new/bias_vizualization
```