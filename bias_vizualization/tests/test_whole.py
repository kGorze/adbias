from pathlib import Path

from bias_vizualization import (
    discover_systems,
    generate_for_systems,
    parse_bias_file,
    renderer_path,
)

EXAMPLE_RESULTS = Path(__file__).parents[1] / "example_results"


def test_whole(tmp_path: Path) -> None:
    systems = discover_systems(EXAMPLE_RESULTS)

    assert len(systems) == 1
    system = systems[0]
    assert system.name == "3CS9"

    biases = parse_bias_file(system.bias_file)
    assert len(biases) == 1
    assert biases[0].bias_type == "don"

    generated = generate_for_systems(
        systems,
        renderer_path(),
        output_directory=tmp_path,
    )

    assert len(generated) == 1
    visualization = generated[0]
    assert visualization.system == "3CS9"
    assert visualization.bias_number == 1
    assert visualization.bias_type == "don"

    output_tcl = tmp_path / "3CS9" / "bias_001_don.tcl"
    assert visualization.output_tcl == output_tcl
    assert output_tcl.is_file()
    assert output_tcl.with_name("bias_001_don_bias.tcl").is_file()
    assert output_tcl.with_name("bias_001_don_scene.pdb").is_file()
