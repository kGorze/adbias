# Przykłady użycia API

## Instalacja

Pakiet wymaga Pythona 3.12 lub nowszego. Z katalogu repozytorium:

```bash
python -m pip install ./benchmark/versions/new/bias_vizualization
```

Podczas rozwijania implementacji warto użyć instalacji edytowalnej. Zmiany w
kodzie będą wtedy widoczne bez ponownej instalacji:

```bash
python -m pip install --editable ./benchmark/versions/new/bias_vizualization
```

Nazwa dystrybucji używana przez `pip` to `bias-vizualization`, a nazwa modułu
w Pythonie to `bias_vizualization`.

## Jedna wizualizacja

`renderer_path()` zwraca zasób Tcl dołączony do paczki, więc aplikacja nie musi
znać położenia kodu źródłowego.

```python
from pathlib import Path

from bias_vizualization import (
    DrawOptions,
    generate_bias_visualization,
    parse_bias_file,
    renderer_path,
)

system = Path("results/3CS9")
bias = parse_bias_file(system / "bias.bpf")[0]

output = generate_bias_visualization(
    mapfile_path=system / "receptor.A.map",
    bias=bias,
    receptor_pdb=system / "receptor_prepared.pdb",
    output_tcl=Path("visualizations/3CS9/bias_001.tcl"),
    renderer_tcl_path=renderer_path(),
    scene_name="3CS9_bias_001",
    epsilon=0.01,
    draw_options=DrawOptions(
        draw_candidate_points=True,
        draw_rejected_points=False,
        accepted_point_radius=0.045,
        graphics_opacity=0.35,
    ),
)
print(output)
```

Powstaną trzy pliki: `bias_001.tcl`, `bias_001_bias.tcl` i
`bias_001_scene.pdb`.

## Jeden jawnie zdefiniowany system

Ten wariant nie zakłada konkretnego układu katalogów wejściowych.

```python
from pathlib import Path

from bias_vizualization import SystemFiles, generate_for_system, renderer_path

system = SystemFiles(
    name="3CS9",
    mapfile=Path("input/3CS9/grid.map"),
    bias_file=Path("input/3CS9/biases.bpf"),
    receptor=Path("input/3CS9/receptor.pdb"),
)

generated = generate_for_system(
    system,
    renderer_path(),
    output_directory=Path("visualizations/3CS9"),
)

for item in generated:
    print(item.bias_number, item.bias_type, item.output_tcl)
```

## Wiele systemów wykrytych w katalogu

`discover_systems()` domyślnie oczekuje plików `receptor.A.map`, `bias.bpf` i
`receptor_prepared.pdb` w osobnym podkatalogu każdego systemu.

```python
from pathlib import Path

from bias_vizualization import (
    discover_systems,
    generate_for_systems,
    renderer_path,
)

systems = discover_systems(Path("results"))
generated = generate_for_systems(
    systems,
    renderer_path(),
    output_directory=Path("visualizations"),
)

for item in generated:
    print(item.system, item.bias_number, item.bias_type, item.output_tcl)
```

Nazwy plików można podać jawnie:

```python
systems = discover_systems(
    "results",
    map_filename="grid.map",
    bias_filename="biases.bpf",
    receptor_filename="protein.pdb",
)
```

## Niższy poziom API

Obliczenia i renderowanie można rozdzielić, aby eksperymentować z budową sceny
bez uruchamiania CLI:

```python
from bias_vizualization import (
    AutoDockGrid,
    DrawOptions,
    build_bias_scene,
    calculate_bias_geometry,
    parse_autodock_mapfile,
    parse_bias_file,
    render_visualization_pdb,
)

spacing, nelements, center = parse_autodock_mapfile("grid.map")
grid = AutoDockGrid(spacing, nelements, center)
bias = parse_bias_file("biases.bpf")[0]
geometry = calculate_bias_geometry(grid, bias, epsilon=0.01)
scene = build_bias_scene(geometry, DrawOptions())
pdb_text = render_visualization_pdb(scene)
```

Stabilnym punktem wejścia dla aplikacji jest `bias_vizualization.api` oraz
symbole eksportowane z `bias_vizualization`. CLI korzysta z tej samej warstwy.
Moduły `generation`, `geometry`, `scene` i `vmd` są implementacją, którą można
zmieniać bez zmiany sygnatur publicznej fasady.

## CLI

Po instalacji dostępna jest komenda:

```bash
bias-vizualization one \
  --map results/3CS9/receptor.A.map \
  --bias-file results/3CS9/bias.bpf \
  --bias-number 1 \
  --receptor results/3CS9/receptor_prepared.pdb \
  --output visualizations/3CS9/bias_001.tcl
```

Dla wszystkich systemów:

```bash
bias-vizualization all \
  --results-dir results \
  --output-directory-name bias_visualizations
```
