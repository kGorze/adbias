from dataclasses import dataclass
import math
from os import PathLike
from pathlib import Path

from .models import Line, Point, Point3D, Primitive, Scene, Sphere


def _format_number(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("cannot write a non-finite number to Tcl")
    return format(value, ".12g")


def _tcl_string(value: str | Path) -> str:
    text = str(value)
    if "\n" in text or "\r" in text:
        raise ValueError("Tcl strings cannot contain line breaks")
    escaped = (
        text.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("$", "\\$")
        .replace("[", "\\[")
        .replace("]", "\\]")
    )
    return f'"{escaped}"'


@dataclass(frozen=True, slots=True)
class _Representation:
    name: str
    label: str
    pdb_resname: str
    color: str
    style: str
    objects: tuple[Primitive, ...]


_GROUP_METADATA = {
    "bias_center": ("Bias center", "BCT"),
    "nearest_grid_point": ("Nearest grid point", "NGP"),
    "center_to_nearest_grid_point": (
        "Bias center to nearest grid point",
        "COF",
    ),
    "one_over_e_bias_surface": ("Bias 1/e isosurface", "BRS"),
    "bias_radius": ("Bias radius", "BRA"),
    "epsilon_energy_surface": ("Epsilon energy isosurface", "EPS"),
    "epsilon_radius": ("Epsilon radius", "EPA"),
    "grid_spacing": ("One grid-spacing step", "GSP"),
    "fixed_candidate_cube": ("Fixed candidate cube", "CUB"),
    "epsilon_grid_box": ("Epsilon-derived grid box", "BOX"),
    "accepted_low_fraction": ("Accepted points: low fraction", "ACY"),
    "accepted_medium_fraction": ("Accepted points: medium fraction", "ACO"),
    "accepted_high_fraction": ("Accepted points: high fraction", "ACR"),
    "rejected_candidate_points": ("Rejected candidate points", "REJ"),
}


def _representations(scene: Scene) -> tuple[_Representation, ...]:
    grouped: dict[str, list[Primitive]] = {}
    for primitive in scene.objects:
        grouped.setdefault(primitive.group, []).append(primitive)

    representations: list[_Representation] = []
    for group, objects in grouped.items():
        try:
            label, pdb_resname = _GROUP_METADATA[group]
        except KeyError as error:
            raise ValueError(f"unknown representation group: {group}") from error

        first = objects[0]
        if any(type(item) is not type(first) for item in objects):
            raise ValueError(f"representation {group} mixes primitive types")
        if any(item.color != first.color for item in objects):
            raise ValueError(f"representation {group} mixes colors")

        match first:
            case Sphere():
                resolution = max(item.resolution for item in objects)
                style = f"VDW 1.0 {resolution}"
            case Line():
                radius = 0.025 * max(item.width for item in objects)
                style = f"Bonds {radius:.3f} 12.0"
            case Point():
                style = "Points 4.0"
            case _:
                raise TypeError(f"unsupported primitive: {type(first).__name__}")

        representations.append(
            _Representation(
                group, label, pdb_resname, first.color, style, tuple(objects)
            )
        )
    return tuple(representations)


def _pdb_atom_line(
    serial: int,
    resname: str,
    resid: int,
    position: Point3D,
    radius: float,
) -> str:
    x, y, z = position
    line = (
        f"HETATM{serial:5d}  V   {resname:>3s} Z{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{radius:6.2f}      "
        "BVIZ V  "
    )
    if len(line) != 80:
        raise ValueError(f"visualization atom does not fit the PDB format: {line}")
    return line


def _line_segments(line: Line) -> tuple[tuple[Point3D, Point3D], ...]:
    if line.style == "solid":
        return ((line.start, line.end),)
    if line.style != "dashed":
        raise ValueError(f"unsupported line style: {line.style}")

    def interpolate(fraction: float) -> Point3D:
        return tuple(
            start + fraction * (end - start)
            for start, end in zip(line.start, line.end, strict=True)
        )

    return tuple(
        (interpolate(index / 16.0), interpolate((index + 1) / 16.0))
        for index in range(0, 16, 2)
    )


def render_visualization_pdb(scene: Scene) -> str:
    """Serialize only the reusable bias geometry, without receptor atoms."""
    serial = 0
    lines: list[str] = []
    bonds: list[tuple[int, int]] = []

    for resid, representation in enumerate(_representations(scene), start=1):
        for primitive in representation.objects:
            match primitive:
                case Sphere(center, radius, _, _, _):
                    serial += 1
                    lines.append(
                        _pdb_atom_line(
                            serial, representation.pdb_resname, resid, center, radius
                        )
                    )
                case Point(position, _, _):
                    serial += 1
                    lines.append(
                        _pdb_atom_line(
                            serial, representation.pdb_resname, resid, position, 0.05
                        )
                    )
                case Line():
                    for start, end in _line_segments(primitive):
                        start_serial = serial + 1
                        end_serial = serial + 2
                        lines.extend(
                            (
                                _pdb_atom_line(
                                    start_serial,
                                    representation.pdb_resname,
                                    resid,
                                    start,
                                    0.05,
                                ),
                                _pdb_atom_line(
                                    end_serial,
                                    representation.pdb_resname,
                                    resid,
                                    end,
                                    0.05,
                                ),
                            )
                        )
                        bonds.append((start_serial, end_serial))
                        serial = end_serial
                case _:
                    raise TypeError(
                        f"unsupported primitive: {type(primitive).__name__}"
                    )
            if serial > 99_999:
                raise ValueError("visualization exceeds the PDB atom serial limit")

    lines.extend(f"CONECT{start:5d}{end:5d}" for start, end in bonds)
    lines.append("END")
    return "\n".join(lines) + "\n"


def render_bias_tcl(
    scene: Scene,
    visualization_pdb: str | PathLike[str],
    renderer_tcl_path: str | PathLike[str],
    bias_name: str,
    graphics_opacity: float,
) -> str:
    """Create a standalone bias loader for any open VMD system."""
    visualization_path = Path(visualization_pdb).resolve()
    renderer_path = Path(renderer_tcl_path).resolve()
    if not visualization_path.is_file():
        raise FileNotFoundError(
            f"visualization PDB does not exist: {visualization_path}"
        )
    if not renderer_path.is_file():
        raise FileNotFoundError(f"renderer Tcl does not exist: {renderer_path}")
    if not math.isfinite(graphics_opacity) or not 0.0 < graphics_opacity <= 1.0:
        raise ValueError("graphics_opacity must be in the range (0, 1]")

    lines = [
        "# Reusable bias overlay. Source this file after loading any system.",
        f"source {_tcl_string(renderer_path)}",
        f"set bias_file {_tcl_string(visualization_path)}",
        (
            "set bias_molecule [::bias_vizualization::load_bias "
            f"$bias_file {_tcl_string(bias_name)} "
            f"{_format_number(graphics_opacity)}]"
        ),
        "",
    ]
    for representation in _representations(scene):
        lines.append(f"# Representation: {representation.label}")
        lines.append(
            f"::bias_vizualization::add_representation $bias_molecule "
            f"{_tcl_string(representation.label)} {representation.pdb_resname} "
            f"{representation.name} {_tcl_string(representation.style)} "
            f"{representation.color}"
        )
    lines.extend(("", "mol top $bias_molecule", ""))
    return "\n".join(lines)


def render_tcl(
    receptor_pdb: str | PathLike[str],
    bias_tcl: str | PathLike[str],
    renderer_tcl_path: str | PathLike[str],
    scene_name: str,
) -> str:
    """Load receptor and source an independently reusable bias overlay."""
    receptor_path = Path(receptor_pdb).resolve()
    bias_tcl_path = Path(bias_tcl).resolve()
    renderer_path = Path(renderer_tcl_path).resolve()
    if not receptor_path.is_file():
        raise FileNotFoundError(f"receptor file does not exist: {receptor_path}")
    if not bias_tcl_path.is_file():
        raise FileNotFoundError(f"bias Tcl does not exist: {bias_tcl_path}")
    if not renderer_path.is_file():
        raise FileNotFoundError(f"renderer Tcl does not exist: {renderer_path}")

    receptor_name = f"{scene_name}_receptor"
    lines = [
        "# Receptor and bias are loaded as separate VMD molecules.",
        f"source {_tcl_string(renderer_path)}",
        f"set receptor_file {_tcl_string(receptor_path)}",
        (
            "set receptor_molecule [::bias_vizualization::load_receptor "
            f"$receptor_file {_tcl_string(receptor_name)}]"
        ),
        f"source {_tcl_string(bias_tcl_path)}",
        "mol top $receptor_molecule",
        "",
    ]
    return "\n".join(lines)
