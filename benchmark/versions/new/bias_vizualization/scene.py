import math
from os import PathLike
from pathlib import Path

from .geometry import add3
from .models import BiasGeometry, DrawOptions, Line, Point, Primitive, Scene, Sphere


def fraction_color(fraction: float) -> str:
    if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be finite and in the range [0, 1]")
    if fraction >= 2.0 / 3.0:
        return "red"
    if fraction >= 1.0 / 3.0:
        return "orange"
    return "yellow"


def build_bias_scene(
    geometry: BiasGeometry,
    mapfile_path: str | PathLike[str],
    draw_options: DrawOptions,
) -> Scene:
    """Convert calculated bias geometry into backend-independent primitives."""
    bias = geometry.bias
    grid = geometry.grid
    center = bias.center
    nearest = geometry.nearest_point
    objects: list[Primitive] = [
        Sphere(center, 0.1, "purple", 20, "bias_center"),
        Sphere(nearest, 0.08, "cyan", 16, "nearest_grid_point"),
    ]

    if geometry.center_grid_distance > 1.0e-9:
        objects.append(
            Line(
                center,
                nearest,
                "cyan",
                3,
                group="center_to_nearest_grid_point",
            )
        )

    radius_end = add3(center, (0.0, bias.radius, 0.0))
    objects.extend(
        (
            Sphere(
                center,
                bias.radius,
                "green",
                30,
                "one_over_e_bias_surface",
            ),
            Line(center, radius_end, "green", 3, group="bias_radius"),
        )
    )

    epsilon_end = add3(center, (0.0, 0.0, geometry.epsilon_radius))
    objects.extend(
        (
            Sphere(
                center,
                geometry.epsilon_radius,
                "red",
                30,
                "epsilon_energy_surface",
            ),
            Line(center, epsilon_end, "red", 3, group="epsilon_radius"),
        )
    )

    point_counts = grid.point_counts
    step_index = (
        geometry.nearest_indices[0] + 1
        if geometry.nearest_indices[0] < point_counts[0] - 1
        else geometry.nearest_indices[0] - 1
    )
    spacing_end = (
        grid.minimum[0] + step_index * grid.spacing,
        nearest[1],
        nearest[2],
    )
    objects.append(Line(nearest, spacing_end, "white", 3, group="grid_spacing"))

    if draw_options.draw_current_cube:
        objects.extend(
            Line(start, end, "gray", group="fixed_candidate_cube")
            for start, end in geometry.current_box_edges
        )
        half_side_end = add3(nearest, (geometry.current_half_side, 0.0, 0.0))
        objects.append(
            Line(nearest, half_side_end, "gray", 3, group="fixed_candidate_cube")
        )

    if draw_options.draw_corrected_box:
        objects.extend(
            Line(start, end, "blue", style="dashed", group="epsilon_grid_box")
            for start, end in geometry.corrected_box_edges
        )

    accepted_point_groups = {
        "yellow": "accepted_low_fraction",
        "orange": "accepted_medium_fraction",
        "red": "accepted_high_fraction",
    }
    for sampled_point in geometry.candidate_points:
        if sampled_point.accepted and draw_options.draw_candidate_points:
            color = fraction_color(sampled_point.fraction)
            objects.append(
                Sphere(
                    sampled_point.position,
                    draw_options.accepted_point_radius,
                    color,
                    4,
                    accepted_point_groups[color],
                )
            )
        elif not sampled_point.accepted and draw_options.draw_rejected_points:
            objects.append(
                Point(sampled_point.position, "gray", "rejected_candidate_points")
            )

    center_status = "yes" if geometry.center_is_on_grid else "no"
    exact_energy_status = (
        "yes"
        if math.isclose(
            geometry.energy_at_nearest_point,
            bias.vset,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        else "no"
    )
    report = (
        "",
        "=== Gaussian bias geometry ===",
        f"map file: {Path(mapfile_path).resolve()}",
        f"grid points: {point_counts[0]} x {point_counts[1]} x {point_counts[2]}",
        f"spacing h: {grid.spacing:.6f} A",
        f"bias center C: {center[0]:.6f} {center[1]:.6f} {center[2]:.6f}",
        (
            f"nearest point G: {nearest[0]:.6f} {nearest[1]:.6f} "
            f"{nearest[2]:.6f}; indices {geometry.nearest_indices[0]} "
            f"{geometry.nearest_indices[1]} {geometry.nearest_indices[2]}"
        ),
        f"distance delta(C,G): {geometry.center_grid_distance:.6f} A",
        f"bias center is on grid: {center_status}",
        (
            f"energy at nearest grid point: "
            f"{geometry.energy_at_nearest_point:.6f} kcal/mol"
        ),
        f"energy at nearest point equals Vset: {exact_energy_status}",
        f"Vset: {bias.vset:.6f} kcal/mol; r: {bias.radius:.6f} A",
        f"fraction at r: exp(-1) = {math.exp(-1.0):.9f}",
        (
            f"epsilon: {geometry.epsilon:.6f} kcal/mol; "
            f"f_epsilon: {geometry.epsilon / abs(bias.vset):.9f}"
        ),
        f"R_epsilon: {geometry.epsilon_radius:.6f} A",
        (
            f"current cube: N={geometry.current_half_intervals}; "
            f"L=N*h={geometry.current_half_side:.6f} A; "
            f"{2 * geometry.current_half_intervals + 1} "
            "points per axis before clipping"
        ),
        f"current candidate points after map clipping: {len(geometry.candidate_points)}",
        f"accepted current-code points: {geometry.accepted_count}",
        (
            f"largest accepted sampled distance: "
            f"{geometry.largest_accepted_distance:.6f} A"
        ),
        "Green sphere: d=r, 1/e isosurface, not a cutoff.",
        "Red sphere: continuous |dE|=epsilon isosurface.",
        "Gray solid box: current fixed 2r-based candidate cube centered at G.",
        "Blue dashed box: grid-aligned bounding box derived from R_epsilon.",
        "",
    )
    return Scene(tuple(objects), report)
