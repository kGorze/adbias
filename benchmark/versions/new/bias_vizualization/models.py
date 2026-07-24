from dataclasses import dataclass
import math
from pathlib import Path


type GridShape = tuple[int, int, int]
type Point3D = tuple[float, float, float]
type AutoDockMapHeader = tuple[float, GridShape, Point3D]
type Index3D = tuple[int, int, int]
type Edge3D = tuple[Point3D, Point3D]


@dataclass(frozen=True, slots=True)
class DrawOptions:
    draw_candidate_points: bool = True
    draw_rejected_points: bool = False
    draw_current_cube: bool = True
    draw_corrected_box: bool = True
    accepted_point_radius: float = 0.045
    graphics_opacity: float = 0.35

    def __post_init__(self) -> None:
        if not math.isfinite(self.accepted_point_radius) or self.accepted_point_radius <= 0.0:
            raise ValueError("accepted_point_radius must be a finite positive number")
        if not math.isfinite(self.graphics_opacity) or not 0.0 < self.graphics_opacity <= 1.0:
            raise ValueError("graphics_opacity must be in the range (0, 1]")


@dataclass(frozen=True, slots=True)
class Bias:
    center: Point3D
    vset: float
    radius: float
    bias_type: str
    source_line: int

    def __post_init__(self) -> None:
        if not all(math.isfinite(value) for value in self.center):
            raise ValueError("bias center coordinates must be finite")
        if not math.isfinite(self.vset) or self.vset >= 0.0:
            raise ValueError("Vset must be a finite negative number")
        if not math.isfinite(self.radius) or self.radius <= 0.0:
            raise ValueError("bias radius must be a finite positive number")
        if not self.bias_type or any(character.isspace() for character in self.bias_type):
            raise ValueError("bias type must be one non-empty token")


@dataclass(frozen=True, slots=True)
class AutoDockGrid:
    spacing: float
    nelements: GridShape
    center: Point3D

    def __post_init__(self) -> None:
        if not math.isfinite(self.spacing) or self.spacing <= 0.0:
            raise ValueError("grid spacing must be a finite positive number")
        if any(value <= 0 for value in self.nelements):
            raise ValueError("grid element counts must be positive")
        if not all(math.isfinite(value) for value in self.center):
            raise ValueError("grid center coordinates must be finite")

    @property
    def point_counts(self) -> GridShape:
        x, y, z = self.nelements
        return x + 1, y + 1, z + 1

    @property
    def minimum(self) -> Point3D:
        cx, cy, cz = self.center
        nx, ny, nz = self.nelements
        half_spacing = self.spacing / 2.0

        return (
            cx - nx * half_spacing,
            cy - ny * half_spacing,
            cz - nz * half_spacing,
        )


@dataclass(frozen=True, slots=True)
class SampledPoint:
    position: Point3D
    distance: float
    fraction: float
    energy_delta: float
    accepted: bool


@dataclass(frozen=True, slots=True)
class BiasGeometry:
    grid: AutoDockGrid
    bias: Bias
    epsilon: float
    nearest_indices: Index3D
    nearest_point: Point3D
    center_grid_distance: float
    energy_at_nearest_point: float
    epsilon_radius: float
    current_half_intervals: int
    current_half_side: float
    current_box_edges: tuple[Edge3D, ...]
    corrected_box_edges: tuple[Edge3D, ...]
    candidate_points: tuple[SampledPoint, ...]
    accepted_count: int
    largest_accepted_distance: float

    @property
    def center_is_on_grid(self) -> bool:
        return self.center_grid_distance <= 1.0e-9


@dataclass(frozen=True, slots=True)
class Sphere:
    center: Point3D
    radius: float
    color: str
    resolution: int = 20


@dataclass(frozen=True, slots=True)
class Line:
    start: Point3D
    end: Point3D
    color: str
    width: int = 2
    style: str = "solid"


@dataclass(frozen=True, slots=True)
class Point:
    position: Point3D
    color: str


@dataclass(frozen=True, slots=True)
class Text:
    position: Point3D
    label: str
    size: float = 1.0


type Primitive = Sphere | Line | Point | Text


@dataclass(frozen=True, slots=True)
class Scene:
    objects: tuple[Primitive, ...]
    report: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GeneratedVisualization:
    system: str
    bias_number: int
    bias_type: str
    output_tcl: Path
