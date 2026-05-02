from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_12_exact
from mms_viscoelastic.visualization import plot_case12_panel, plot_scalar_field, save_case12_component_maps


def test_plot_scalar_field_saves_file(tmp_path: Path) -> None:
    grid = Grid2D.from_bounds(nx=9, ny=7)
    field = evaluate_case_12_exact(grid)
    output = tmp_path / "u.png"

    fig, ax = plot_scalar_field(grid, field.u, title="u", output_path=output)
    plt.close(fig)

    assert output.exists()
    assert output.stat().st_size > 0
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"


def test_plot_scalar_field_rejects_swapped_axes() -> None:
    grid = Grid2D.from_bounds(nx=9, ny=7)
    field = evaluate_case_12_exact(grid)

    with pytest.raises(ValueError, match="shape"):
        plot_scalar_field(grid, field.u.T, title="u transposed")


def test_case12_panel_and_component_maps(tmp_path: Path) -> None:
    grid = Grid2D.from_bounds(nx=11, ny=11)
    field = evaluate_case_12_exact(grid)

    panel_path = tmp_path / "panel.png"
    panel = plot_case12_panel(grid, field, output_path=panel_path)
    plt.close(panel)

    paths = save_case12_component_maps(grid, field, tmp_path / "components", components=["u", "v"])

    assert panel_path.exists()
    assert len(paths) == 2
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)
