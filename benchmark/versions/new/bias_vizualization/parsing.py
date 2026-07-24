import math
from os import PathLike
from pathlib import Path

from .models import AutoDockMapHeader, Bias, GridShape, Point3D


def parse_autodock_mapfile(
    mapfile_path: str | PathLike[str],
) -> AutoDockMapHeader:
    """Read SPACING, NELEMENTS and CENTER from the first six map lines."""
    path = Path(mapfile_path)
    spacing: float | None = None
    nelements: GridShape | None = None
    center: Point3D | None = None

    with path.open("r", encoding="utf-8") as file:
        for line_number in range(1, 7):
            line = file.readline()
            if line == "":
                raise ValueError(
                    f"{path}: expected six header lines, found {line_number - 1}"
                )

            fields = line.split()
            if not fields or fields[0] not in {"SPACING", "NELEMENTS", "CENTER"}:
                continue

            key = fields[0]
            expected_values = 1 if key == "SPACING" else 3
            if len(fields) != expected_values + 1:
                raise ValueError(
                    f"{path}:{line_number}: {key} requires exactly "
                    f"{expected_values} value(s), found {len(fields) - 1}"
                )

            values = fields[1:]
            try:
                if key == "SPACING":
                    if spacing is not None:
                        raise ValueError("duplicate SPACING field")
                    spacing = float(values[0])
                    if not math.isfinite(spacing) or spacing <= 0.0:
                        raise ValueError("SPACING must be a finite positive number")
                elif key == "NELEMENTS":
                    if nelements is not None:
                        raise ValueError("duplicate NELEMENTS field")
                    nelements = int(values[0]), int(values[1]), int(values[2])
                    if any(value <= 0 for value in nelements):
                        raise ValueError("NELEMENTS values must be positive integers")
                else:
                    if center is not None:
                        raise ValueError("duplicate CENTER field")
                    center = float(values[0]), float(values[1]), float(values[2])
                    if not all(math.isfinite(value) for value in center):
                        raise ValueError("CENTER values must be finite numbers")
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid {key}: {error}"
                ) from error

    missing: list[str] = []
    if spacing is None:
        missing.append("SPACING")
    if nelements is None:
        missing.append("NELEMENTS")
    if center is None:
        missing.append("CENTER")
    if missing:
        raise ValueError(
            f"{path}: missing {', '.join(missing)} in the first six header lines"
        )

    return spacing, nelements, center


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
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected 6 fields, found {len(fields)}"
                )
            try:
                bias = Bias(
                    center=(float(fields[0]), float(fields[1]), float(fields[2])),
                    vset=float(fields[3]),
                    radius=float(fields[4]),
                    bias_type=fields[5],
                    source_line=line_number,
                )
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid bias: {error}") from error
            biases.append(bias)

    if not biases:
        raise ValueError(f"{path}: no bias definitions found")
    return tuple(biases)
