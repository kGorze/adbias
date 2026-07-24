from collections.abc import Sequence
from os import PathLike
from pathlib import Path

from .geometry import calculate_bias_geometry
from .models import AutoDockGrid, Bias, DrawOptions, GeneratedVisualization
from .parsing import parse_autodock_mapfile, parse_bias_file
from .scene import build_bias_scene
from .vmd import render_bias_tcl, render_tcl, render_visualization_pdb


def generate_bias_visualization(
    mapfile_path: str | PathLike[str],
    bias: Bias,
    receptor_pdb: str | PathLike[str],
    output_tcl: str | PathLike[str],
    renderer_tcl_path: str | PathLike[str],
    scene_name: str,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
) -> Path:
    spacing, nelements, center = parse_autodock_mapfile(mapfile_path)
    grid = AutoDockGrid(spacing, nelements, center)
    geometry = calculate_bias_geometry(grid, bias, epsilon)

    scene = build_bias_scene(geometry, mapfile_path, draw_options)

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
        raise ValueError("bias type does not contain filename-safe characters")
    return token


def generate_for_system(
    system_dir: str | PathLike[str],
    renderer_tcl_path: str | PathLike[str],
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory_name: str = "bias_visualizations",
    map_filename: str = "receptor.A.map",
    bias_filename: str = "bias.bpf",
    receptor_filename: str = "receptor_prepared.pdb",
) -> tuple[GeneratedVisualization, ...]:
    directory = Path(system_dir)
    if not directory.is_dir():
        raise FileNotFoundError(f"system directory does not exist: {directory}")

    mapfile_path = directory / map_filename
    bias_file_path = directory / bias_filename
    receptor_path = directory / receptor_filename
    for required_path in (mapfile_path, bias_file_path, receptor_path):
        if not required_path.is_file():
            raise FileNotFoundError(f"required system file does not exist: {required_path}")

    biases = parse_bias_file(bias_file_path)
    output_directory = directory / output_directory_name
    generated: list[GeneratedVisualization] = []
    for bias_number, bias in enumerate(biases, start=1):
        bias_type = _safe_filename_token(bias.bias_type)
        output_tcl = output_directory / f"bias_{bias_number:03d}_{bias_type}.tcl"
        scene_name = f"{directory.name}_bias_{bias_number:03d}_{bias_type}"
        generate_bias_visualization(
            mapfile_path=mapfile_path,
            bias=bias,
            receptor_pdb=receptor_path,
            output_tcl=output_tcl,
            renderer_tcl_path=renderer_tcl_path,
            scene_name=scene_name,
            epsilon=epsilon,
            draw_options=draw_options,
        )
        generated.append(
            GeneratedVisualization(
                system=directory.name,
                bias_number=bias_number,
                bias_type=bias.bias_type,
                output_tcl=output_tcl,
            )
        )
    return tuple(generated)


def generate_for_all_systems(
    results_directory: str | PathLike[str],
    renderer_tcl_path: str | PathLike[str],
    systems: Sequence[str] | None = None,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory_name: str = "bias_visualizations",
    map_filename: str = "receptor.A.map",
    bias_filename: str = "bias.bpf",
    receptor_filename: str = "receptor_prepared.pdb",
) -> tuple[GeneratedVisualization, ...]:
    results_path = Path(results_directory)
    if not results_path.is_dir():
        raise FileNotFoundError(f"results directory does not exist: {results_path}")

    if systems is None:
        system_names = tuple(
            directory.name
            for directory in sorted(results_path.iterdir())
            if directory.is_dir() and (directory / bias_filename).is_file()
        )
        if not system_names:
            raise ValueError(f"no system directories with {bias_filename} in {results_path}")
    else:
        system_names = tuple(systems)
        if not system_names:
            raise ValueError("systems cannot be empty")

    generated: list[GeneratedVisualization] = []
    for system_name in system_names:
        generated.extend(
            generate_for_system(
                system_dir=results_path / system_name,
                renderer_tcl_path=renderer_tcl_path,
                epsilon=epsilon,
                draw_options=draw_options,
                output_directory_name=output_directory_name,
                map_filename=map_filename,
                bias_filename=bias_filename,
                receptor_filename=receptor_filename,
            )
        )
    return tuple(generated)
