from collections.abc import Sequence
from pathlib import Path

from .generation import (
    SystemFiles,
    discover_systems as _discover_systems,
    generate_bias_visualization as _generate_bias_visualization,
    generate_for_system as _generate_for_system,
    generate_for_systems as _generate_for_systems,
)
from .models import Bias, DrawOptions, GeneratedVisualization
from .parsing import parse_bias_file as _parse_bias_file


def renderer_path() -> Path:
    """Return the VMD renderer installed with the package."""
    path = Path(__file__).with_name("renderer.tcl")
    if not path.is_file():
        raise FileNotFoundError(f"package renderer does not exist: {path}")
    return path


def parse_bias_file(bias_file_path: str | Path) -> tuple[Bias, ...]:
    """Read all bias definitions from a BPF file."""
    return _parse_bias_file(bias_file_path)


def discover_systems(
    results_directory: str | Path,
    *,
    map_filename: str = "receptor.A.map",
    bias_filename: str = "bias.bpf",
    receptor_filename: str = "receptor_prepared.pdb",
) -> tuple[SystemFiles, ...]:
    """Discover systems using the standard results-directory layout."""
    return _discover_systems(
        results_directory,
        map_filename=map_filename,
        bias_filename=bias_filename,
        receptor_filename=receptor_filename,
    )


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
    """Generate one complete VMD visualization."""
    return _generate_bias_visualization(
        mapfile_path=mapfile_path,
        bias=bias,
        receptor_pdb=receptor_pdb,
        output_tcl=output_tcl,
        renderer_tcl_path=renderer_tcl_path,
        scene_name=scene_name,
        epsilon=epsilon,
        draw_options=draw_options,
    )


def generate_for_system(
    system: SystemFiles,
    renderer_tcl_path: str | Path,
    *,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory: str | Path,
) -> tuple[GeneratedVisualization, ...]:
    """Generate visualizations for every bias in one system."""
    return _generate_for_system(
        system,
        renderer_tcl_path,
        epsilon=epsilon,
        draw_options=draw_options,
        output_directory=output_directory,
    )


def generate_for_systems(
    systems: Sequence[SystemFiles],
    renderer_tcl_path: str | Path,
    *,
    epsilon: float = 0.01,
    draw_options: DrawOptions = DrawOptions(),
    output_directory: str | Path,
) -> tuple[GeneratedVisualization, ...]:
    """Generate visualizations for every bias in multiple systems."""
    return _generate_for_systems(
        systems,
        renderer_tcl_path,
        epsilon=epsilon,
        draw_options=draw_options,
        output_directory=output_directory,
    )
