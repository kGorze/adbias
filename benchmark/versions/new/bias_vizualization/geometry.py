import math

from .models import AutoDockGrid, Bias, BiasGeometry, Edge3D, Point3D, SampledPoint


def _clamp_int(value: int, minimum: int, maximum: int) -> int:
    if minimum > maximum:
        raise ValueError("minimum cannot be greater than maximum")
    return max(minimum, min(value, maximum))


def nearest_grid_index(
    coordinate: float,
    minimum: float,
    spacing: float,
    point_count: int,
) -> int:
    if spacing <= 0.0:
        raise ValueError("spacing must be positive")
    if point_count <= 0:
        raise ValueError("point_count must be positive")
    unbounded_index = math.floor((coordinate - minimum) / spacing + 0.5)
    return _clamp_int(unbounded_index, 0, point_count - 1)


def distance3(first: Point3D, second: Point3D) -> float:
    return math.dist(first, second)


def add3(first: Point3D, second: Point3D) -> Point3D:
    return tuple(
        left + right for left, right in zip(first, second, strict=True)
    )


def _box_edges(minimum: Point3D, maximum: Point3D) -> tuple[Edge3D, ...]:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    vertices: tuple[Point3D, ...] = (
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    )
    edge_indices = (
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    )
    return tuple((vertices[start], vertices[end]) for start, end in edge_indices)


def _corrected_axis_range(
    bias_coordinate: float,
    epsilon_radius: float,
    grid_minimum: float,
    spacing: float,
    point_count: int,
) -> tuple[int, int] | None:
    tolerance = 1.0e-12
    first = math.ceil(
        (bias_coordinate - epsilon_radius - grid_minimum) / spacing - tolerance
    )
    last = math.floor(
        (bias_coordinate + epsilon_radius - grid_minimum) / spacing + tolerance
    )
    first = max(first, 0)
    last = min(last, point_count - 1)
    if first > last:
        return None
    return first, last


def calculate_bias_geometry(
    grid: AutoDockGrid,
    bias: Bias,
    epsilon: float = 0.01,
) -> BiasGeometry:
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be a finite positive number")

    amplitude = abs(bias.vset)
    if amplitude <= epsilon:
        raise ValueError("abs(Vset) must be greater than epsilon")

    point_counts = grid.point_counts
    grid_minimum = grid.minimum
    nearest_indices = tuple(
        nearest_grid_index(coordinate, minimum, grid.spacing, point_count)
        for coordinate, minimum, point_count in zip(
            bias.center,
            grid_minimum,
            point_counts,
            strict=True,
        )
    )
    nearest_point = tuple(
        minimum + index * grid.spacing
        for minimum, index in zip(grid_minimum, nearest_indices, strict=True)
    )
    center_grid_distance = distance3(bias.center, nearest_point)
    energy_at_nearest_point = bias.vset * math.exp(
        -(center_grid_distance * center_grid_distance) / (bias.radius * bias.radius)
    )
    epsilon_radius = bias.radius * math.sqrt(math.log(amplitude / epsilon))

    half_intervals = math.ceil(2.0 * bias.radius / grid.spacing)
    half_side = half_intervals * grid.spacing
    current_minimum = tuple(value - half_side for value in nearest_point)
    current_maximum = tuple(value + half_side for value in nearest_point)
    current_box_edges = _box_edges(current_minimum, current_maximum)

    current_ranges = tuple(
        (
            _clamp_int(index - half_intervals, 0, point_count - 1),
            _clamp_int(index + half_intervals, 0, point_count - 1),
        )
        for index, point_count in zip(nearest_indices, point_counts, strict=True)
    )

    corrected_ranges = tuple(
        _corrected_axis_range(
            coordinate,
            epsilon_radius,
            minimum,
            grid.spacing,
            point_count,
        )
        for coordinate, minimum, point_count in zip(
            bias.center,
            grid_minimum,
            point_counts,
            strict=True,
        )
    )
    if any(axis_range is None for axis_range in corrected_ranges):
        corrected_box_edges: tuple[Edge3D, ...] = ()
    else:
        complete_ranges = tuple(
            axis_range for axis_range in corrected_ranges if axis_range is not None
        )
        corrected_minimum = tuple(
            minimum + axis_range[0] * grid.spacing
            for minimum, axis_range in zip(
                grid_minimum,
                complete_ranges,
                strict=True,
            )
        )
        corrected_maximum = tuple(
            minimum + axis_range[1] * grid.spacing
            for minimum, axis_range in zip(
                grid_minimum,
                complete_ranges,
                strict=True,
            )
        )
        corrected_box_edges = _box_edges(corrected_minimum, corrected_maximum)

    x_range, y_range, z_range = current_ranges
    bx, by, bz = bias.center
    candidate_points: list[SampledPoint] = []
    accepted_count = 0
    largest_accepted_distance = 0.0
    for k in range(z_range[0], z_range[1] + 1):
        z = grid_minimum[2] + k * grid.spacing
        for j in range(y_range[0], y_range[1] + 1):
            y = grid_minimum[1] + j * grid.spacing
            for i in range(x_range[0], x_range[1] + 1):
                x = grid_minimum[0] + i * grid.spacing
                squared_distance = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
                distance = math.sqrt(squared_distance)
                fraction = math.exp(-squared_distance / (bias.radius * bias.radius))
                energy_delta = bias.vset * fraction
                accepted = energy_delta < -epsilon
                candidate_points.append(
                    SampledPoint(
                        position=(x, y, z),
                        distance=distance,
                        fraction=fraction,
                        energy_delta=energy_delta,
                        accepted=accepted,
                    )
                )
                if accepted:
                    accepted_count += 1
                    largest_accepted_distance = max(largest_accepted_distance, distance)

    return BiasGeometry(
        grid=grid,
        bias=bias,
        epsilon=epsilon,
        nearest_indices=nearest_indices,
        nearest_point=nearest_point,
        center_grid_distance=center_grid_distance,
        energy_at_nearest_point=energy_at_nearest_point,
        epsilon_radius=epsilon_radius,
        current_half_intervals=half_intervals,
        current_half_side=half_side,
        current_box_edges=current_box_edges,
        corrected_box_edges=corrected_box_edges,
        candidate_points=tuple(candidate_points),
        accepted_count=accepted_count,
        largest_accepted_distance=largest_accepted_distance,
    )
