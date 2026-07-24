import math

from .models import BiasGeometry, DrawOptions, Line, Point, Primitive, Scene, Sphere


#frakcje są arbitralne i można je definiować na wiele sposobów
def fraction_color(fraction: float) -> str:
    if fraction >= 2.0 / 3.0:
        return "red"
    if fraction >= 1.0 / 3.0:
        return "orange"
    return "yellow"


def build_bias_scene(
    geometry: BiasGeometry,
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

    objects.extend(
        (
            Sphere(
                center,
                bias.radius,
                "green",
                30,
                "one_over_e_bias_surface",
            ),
            Sphere(
                center,
                geometry.epsilon_radius,
                "red",
                30,
                "epsilon_energy_surface",
            ),
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

    return Scene(tuple(objects))
