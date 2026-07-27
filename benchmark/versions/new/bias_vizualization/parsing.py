import math
from itertools import islice
from os import PathLike
from pathlib import Path

from .models import AutoDockMapHeader, Bias, GridShape, Point3D



#parser jest dostosowany idealnie pod format pliku mapfile. tutaj jest przykład
"""
GRID_PARAMETER_FILE receptor.gpf
GRID_DATA_FILE receptor.maps.fld
MACROMOLECULE receptor.pdbqt
SPACING 0.375
NELEMENTS 90 72 68
CENTER 29.025 2.847 52.761
0
0
...
"""
#TODO: można zrobić to z większą starannością w sprawdzaniu pliku, żeby format był w stu procentach pewny
def parse_autodock_mapfile(
    mapfile_path: str | PathLike[str],
) -> AutoDockMapHeader:
    """Read SPACING, NELEMENTS and CENTER from fixed map header lines 4-6."""
    path = Path(mapfile_path)

    with path.open("r", encoding="utf-8") as file:
        header_lines = tuple(islice(file, 6))

    if len(header_lines) != 6:
        raise ValueError(f"{path}: expected six header lines, found {len(header_lines)}")

    spacing_fields = header_lines[3].split()
    nelements_fields = header_lines[4].split()
    center_fields = header_lines[5].split()

    if len(spacing_fields) != 2 or spacing_fields[0] != "SPACING":
        raise ValueError(f"{path}:4: expected 'SPACING value'")
    if len(nelements_fields) != 4 or nelements_fields[0] != "NELEMENTS":
        raise ValueError(f"{path}:5: expected 'NELEMENTS x y z'")
    if len(center_fields) != 4 or center_fields[0] != "CENTER":
        raise ValueError(f"{path}:6: expected 'CENTER x y z'")

    spacing = float(spacing_fields[1])
    nelements: GridShape = (
        int(nelements_fields[1]),
        int(nelements_fields[2]),
        int(nelements_fields[3]),
    )
    center: Point3D = (
        float(center_fields[1]),
        float(center_fields[2]),
        float(center_fields[3]),
    )

    #wartości muszą być w tych granicach
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError(f"{path}:4: SPACING must be a finite positive number")
    if any(value <= 0 for value in nelements):
        raise ValueError(f"{path}:5: NELEMENTS values must be positive integers")
    if not all(math.isfinite(value) for value in center):
        raise ValueError(f"{path}:6: CENTER values must be finite numbers")

    return spacing, nelements, center


#również jest bardzo sztywno zrobione, format musi być i powinien być dokładnie taki jaki jest w zazwyczaj
"""
x y z Vset r type
30.620 -0.693 52.109 -1.50 1.00 don
"""
def parse_bias_file(bias_file_path: str | PathLike[str]) -> tuple[Bias, ...]:
    path = Path(bias_file_path)
    with path.open("r", encoding="utf-8") as file:
        header = file.readline().split()
        expected_header = ["x", "y", "z", "Vset", "r", "type"]
        if header != expected_header:
            raise ValueError(
                f"{path}:1: expected header {' '.join(expected_header)!r}, "
                f"found {' '.join(header)!r}"
            )

        biases: list[Bias] = []
        for line_number, line in enumerate(file, start=2):
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected 6 fields, found {len(fields)}"
                )
            biases.append(
                Bias(
                    center=(float(fields[0]), float(fields[1]), float(fields[2])),
                    vset=float(fields[3]),
                    radius=float(fields[4]),
                    bias_type=fields[5],
                    source_line=line_number,
                )
            )

    if not biases:
        raise ValueError(f"{path}: no bias definitions found")
    return tuple(biases)
