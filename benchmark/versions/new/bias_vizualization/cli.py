import argparse
from pathlib import Path

from .generation import generate_bias_visualization, generate_for_all_systems
from .models import DrawOptions
from .parsing import parse_bias_file


def _add_drawing_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--opacity", type=float, default=0.35)
    parser.add_argument("--accepted-point-radius", type=float, default=0.045)
    parser.add_argument("--draw-rejected-points", action="store_true")
    parser.add_argument("--no-candidate-points", action="store_true")
    parser.add_argument("--no-current-cube", action="store_true")
    parser.add_argument("--no-corrected-box", action="store_true")


def _draw_options_from_arguments(arguments: argparse.Namespace) -> DrawOptions:
    return DrawOptions(
        draw_candidate_points=not arguments.no_candidate_points,
        draw_rejected_points=arguments.draw_rejected_points,
        draw_current_cube=not arguments.no_current_cube,
        draw_corrected_box=not arguments.no_corrected_box,
        accepted_point_radius=arguments.accepted_point_radius,
        graphics_opacity=arguments.opacity,
    )


def _argument_parser() -> argparse.ArgumentParser:
    package_directory = Path(__file__).resolve().parent
    default_renderer = package_directory / "renderer.tcl"
    default_results = package_directory.parent / "results"
    parser = argparse.ArgumentParser(
        description="Calculate AutoDock bias geometry and generate drawing-only VMD Tcl files."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    one_parser = subparsers.add_parser("one", help="generate Tcl for one bias")
    one_parser.add_argument("--map", required=True, dest="mapfile", type=Path)
    one_parser.add_argument("--bias-file", required=True, type=Path)
    one_parser.add_argument("--bias-number", type=int, default=1)
    one_parser.add_argument("--receptor", required=True, type=Path)
    one_parser.add_argument("--output", required=True, type=Path)
    one_parser.add_argument("--scene-name", default="Gaussian_bias")
    one_parser.add_argument("--renderer", type=Path, default=default_renderer)
    _add_drawing_arguments(one_parser)

    all_parser = subparsers.add_parser(
        "all",
        help="generate one Tcl file per bias for every requested system",
    )
    all_parser.add_argument("--results-dir", type=Path, default=default_results)
    all_parser.add_argument(
        "--systems",
        nargs="+",
        help="system names; omit to discover systems in the results directory",
    )
    all_parser.add_argument("--renderer", type=Path, default=default_renderer)
    all_parser.add_argument("--output-directory-name", default="bias_visualizations")
    all_parser.add_argument("--map-filename", default="receptor.A.map")
    all_parser.add_argument("--bias-filename", default="bias.bpf")
    all_parser.add_argument("--receptor-filename", default="receptor_prepared.pdb")
    _add_drawing_arguments(all_parser)
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    draw_options = _draw_options_from_arguments(arguments)

    if arguments.command == "one":
        biases = parse_bias_file(arguments.bias_file)
        if not 1 <= arguments.bias_number <= len(biases):
            raise ValueError(
                f"bias number {arguments.bias_number} is outside 1..{len(biases)}"
            )
        output = generate_bias_visualization(
            mapfile_path=arguments.mapfile,
            bias=biases[arguments.bias_number - 1],
            receptor_pdb=arguments.receptor,
            output_tcl=arguments.output,
            renderer_tcl_path=arguments.renderer,
            scene_name=arguments.scene_name,
            epsilon=arguments.epsilon,
            draw_options=draw_options,
        )
        print(output)
        return

    generated = generate_for_all_systems(
        results_directory=arguments.results_dir,
        renderer_tcl_path=arguments.renderer,
        systems=arguments.systems,
        epsilon=arguments.epsilon,
        draw_options=draw_options,
        output_directory_name=arguments.output_directory_name,
        map_filename=arguments.map_filename,
        bias_filename=arguments.bias_filename,
        receptor_filename=arguments.receptor_filename,
    )
    for visualization in generated:
        print(
            f"{visualization.system} bias {visualization.bias_number} "
            f"({visualization.bias_type}): {visualization.output_tcl}"
        )
