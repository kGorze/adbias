from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from .geometry import calculate_bias_geometry
from .models import AutoDockGrid, Bias, DrawOptions, GeneratedVisualization
from .parsing import parse_autodock_mapfile, parse_bias_file
from .scene import build_bias_scene
from .vmd import render_bias_tcl, render_tcl, render_visualization_pdb


@dataclass(frozen=True, slots=True)
class SystemFiles:
    mapfile: Path
    bias_file: Path
    receptor: Path
    name: str


def generate_bias_visualization(
    mapfile_path: str | Path,
    bias: Bias,
    receptor_pdb: str | Path,
    output_tcl: str | Path,
    renderer_tcl_path: str | Path,
    scene_name: str,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
) -> Path:
    spacing, nelements, center = parse_autodock_mapfile(mapfile_path)
    grid = AutoDockGrid(spacing, nelements, center)

    geometry = calculate_bias_geometry(grid, bias, epsilon)
    scene = build_bias_scene(geometry, draw_options)

    output_path = Path(output_tcl)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    visualization_path = output_path.with_name(f"{output_path.stem}_scene.pdb")
    visualization_path.write_text(
        render_visualization_pdb(scene),
        encoding="utf-8",
    )
    bias_tcl_path = output_path.with_name(f"{output_path.stem}_bias.tcl")
    bias_tcl_path.write_text(
        render_bias_tcl(
            scene=scene,
            visualization_pdb=visualization_path,
            renderer_tcl_path=renderer_tcl_path,
            bias_name=f"{scene_name}_bias",
            graphics_opacity=draw_options.graphics_opacity,
        ),
        encoding="utf-8",
    )
    tcl_script = render_tcl(
        receptor_pdb=receptor_pdb,
        bias_tcl=bias_tcl_path,
        renderer_tcl_path=renderer_tcl_path,
        scene_name=scene_name,
    )
    output_path.write_text(tcl_script, encoding="utf-8")
    return output_path


def _safe_filename_token(value: str) -> str:
    token = "".join(
        character if character.isalnum() or character in {"-", "_"} else "_"
        for character in value
    )
    if not token:
        raise ValueError("value does not contain filename-safe characters")
    return token


def generate_for_system(
    system: SystemFiles,
    renderer_tcl_path: str | Path,
    *,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory: str | Path,
) -> tuple[GeneratedVisualization, ...]:
    """Generate visualizations for one explicitly defined system."""
    if not system.name:
        raise ValueError("system name cannot be empty")

    for required_path in (system.mapfile, system.bias_file, system.receptor):
        if not required_path.is_file():
            raise FileNotFoundError(f"required system file does not exist: {required_path}")

    biases = parse_bias_file(system.bias_file)
    output_path = Path(output_directory)
    generated: list[GeneratedVisualization] = []
    for bias_number, bias in enumerate(biases, start=1):
        bias_type = _safe_filename_token(bias.bias_type)
        output_tcl = output_path / f"bias_{bias_number:03d}_{bias_type}.tcl"
        scene_name = f"{system.name}_bias_{bias_number:03d}_{bias_type}"
        generate_bias_visualization(
            mapfile_path=system.mapfile,
            bias=bias,
            receptor_pdb=system.receptor,
            output_tcl=output_tcl,
            renderer_tcl_path=renderer_tcl_path,
            scene_name=scene_name,
            epsilon=epsilon,
            draw_options=draw_options,
        )
        generated.append(
            GeneratedVisualization(
                system=system.name,
                bias_number=bias_number,
                bias_type=bias.bias_type,
                output_tcl=output_tcl,
            )
        )
    return tuple(generated)


def generate_for_systems(
    systems: Sequence[SystemFiles],
    renderer_tcl_path: str | Path,
    *,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory: str | Path,
) -> tuple[GeneratedVisualization, ...]:
    """Generate visualizations in one output subdirectory per system."""
    if not systems:
        raise ValueError("systems cannot be empty")

    output_directories: dict[str, str] = {}
    for system in systems:
        directory_name = _safe_filename_token(system.name)
        conflicting_system = output_directories.get(directory_name)
        if conflicting_system is not None:
            raise ValueError(
                f"systems {conflicting_system!r} and {system.name!r} "
                f"map to the same output directory {directory_name!r}"
            )
        output_directories[directory_name] = system.name

    output_path = Path(output_directory)
    generated: list[GeneratedVisualization] = []
    for system in systems:
        generated.extend(
            generate_for_system(
                system,
                renderer_tcl_path,
                epsilon=epsilon,
                draw_options=draw_options,
                output_directory=output_path / _safe_filename_token(system.name),
            )
        )

    return tuple(generated)


def discover_systems(
    results_directory: str | Path,
    *,
    map_filename: str = "receptor.A.map",
    bias_filename: str = "bias.bpf",
    receptor_filename: str = "receptor_prepared.pdb",
) -> tuple[SystemFiles, ...]:
    """Discover systems using the standard results-directory layout."""
    results_path = Path(results_directory)
    if not results_path.is_dir():
        raise FileNotFoundError(f"results directory does not exist: {results_path}")

    systems = tuple(
        SystemFiles(
            name=directory.name,
            mapfile=directory / map_filename,
            bias_file=directory / bias_filename,
            receptor=directory / receptor_filename,
        )
        for directory in sorted(results_path.iterdir())

        #system dopiero jest wykryty, jeżeli jest plik bias.bpf w katalogu, inne nie są brane pod uwagę
        #TODO: można zmienić, żeby sprawdzało resztę plików. to jest zależne od preferencji, jednak to co robi się z receptorem albo mapą nie jest dla nas istotne
        if directory.is_dir() and (directory / bias_filename).is_file()
    )
    if not systems:
        raise ValueError(
            f"no systems containing {bias_filename!r} found in {results_path}"
        )

    return systems
