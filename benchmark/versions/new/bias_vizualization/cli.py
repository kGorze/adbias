import argparse
from pathlib import Path

from .generation import (
    discover_systems,
    generate_bias_visualization,
    generate_for_systems,
)
from .models import DrawOptions
from .parsing import parse_bias_file


def _add_drawing_arguments(parser: argparse.ArgumentParser) -> None:
    # wspólne argumenty dotyczące wyglądu wizualizacji dla komend "one" i "all"
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--opacity", type=float, default=0.35)
    parser.add_argument("--accepted-point-radius", type=float, default=0.045)
    parser.add_argument("--draw-rejected-points", action="store_true")
    parser.add_argument("--no-candidate-points", action="store_true")


def _draw_options_from_arguments(arguments: argparse.Namespace) -> DrawOptions:
    # zamiana argumentów argparse na jeden obiekt z opcjami rysowania
    # no_candidate_points trzeba odwrócić, ponieważ DrawOptions używa draw_candidate_points
    return DrawOptions(
        draw_candidate_points=not arguments.no_candidate_points,
        draw_rejected_points=arguments.draw_rejected_points,
        accepted_point_radius=arguments.accepted_point_radius,
        graphics_opacity=arguments.opacity,
    )


def _argument_parser() -> argparse.ArgumentParser:
    package_directory = Path(__file__).resolve().parent
    default_renderer = package_directory / "renderer.tcl"
    default_results = package_directory.parent / "results"
    parser = argparse.ArgumentParser()
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


    # nargs="+" oznacza co najmniej jedną nazwę systemu, jeżeli --systems zostanie podane
    # brak --systems oznacza, że użyte zostaną wszystkie wykryte systemy
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
    args = _argument_parser().parse_args()

    # wspólne argumenty rysowania składamy w jeden DrawOptions
    draw_options = _draw_options_from_arguments(args)

    if args.command == "one":
        # jeden plik bias.bpf może zawierać wiele biasów
        #TODO: warto się tym zająć i zrobić testy co w takim wypadku, bo może być problem
        biases = parse_bias_file(args.bias_file)

        if not 1 <= args.bias_number <= len(biases):
            raise ValueError(
                f"bias number {args.bias_number} is outside 1..{len(biases)}"
            )

        output = generate_bias_visualization(
            mapfile_path=args.mapfile,
            # -1 zamienia numer biasu podany przez użytkownika na indeks listy
            bias=biases[args.bias_number - 1],
            receptor_pdb=args.receptor,
            output_tcl=args.output,
            renderer_tcl_path=args.renderer,
            scene_name=args.scene_name,
            epsilon=args.epsilon,
            draw_options=draw_options,
        )

        print(output)
        return

    # discover_systems wykrywa systemy na podstawie konkretnego układu katalogów
    # i nazw plików, więc zmiana struktury results może spowodować brak wykrytych systemów
    systems = discover_systems(
        args.results_dir,
        map_filename=args.map_filename,
        bias_filename=args.bias_filename,
        receptor_filename=args.receptor_filename,
    )

    # jeżeli użytkownik podał --systems, wybieramy tylko wskazane systemy
    if args.systems is not None:
        systems_by_name = {system.name: system for system in systems}
        missing = [name for name in args.systems if name not in systems_by_name]

        # sprawdzamy wcześniej, czy użytkownik nie podał nieistniejącego systemu
        if missing:
            raise ValueError("systems not found: " + ", ".join(missing))

        # zachowujemy kolejność systemów w --systems
        systems = tuple(systems_by_name[name] for name in args.systems)

    # dla każdego systemu i każdego biasu generowany jest osobny plik Tcl
    generated = generate_for_systems(
        systems,
        args.renderer,
        epsilon=args.epsilon,
        draw_options=draw_options,

        # wszystkie wizualizacje trafiaja do jednego katalogu wewnatrz results
        output_directory=args.results_dir / args.output_directory_name,
    )

    # na koncu wypisujemy informacje o kazdym wygenerowanym pliku
    for visualization in generated:
        print(
            f"{visualization.system} bias {visualization.bias_number} "
            f"({visualization.bias_type}): {visualization.output_tcl}"
        )
