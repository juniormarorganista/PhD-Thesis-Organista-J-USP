from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_toy_solution


def test_grid_shape_and_field_features() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_toy_solution(grid)

    assert grid.X.shape == (4, 5)
    assert grid.Y.shape == (4, 5)
    assert field.as_feature_matrix().shape == (20, 7)
