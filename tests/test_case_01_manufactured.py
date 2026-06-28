import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case01Parameters,
    evaluate_case_01_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_01_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_01_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_01_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_01_exact(grid, t=0.0, parameters=Case01Parameters(beta=0.5))

    # Point (x, y) = (0, 0): the velocity/vorticity/stream-function fields
    # vanish at x=0 (factor (x-1)^2 x^2), and the stresses are spatially
    # uniform: (1 - beta) * (c1, c2, c3) = 0.5 * (1, 2, 3).
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 1.0)
    assert np.isclose(field.tyy[0, 0], 1.5)


def test_case_01_stress_field_is_spatially_uniform() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=5)
    field = evaluate_case_01_exact(grid, t=0.0, parameters=Case01Parameters(beta=0.5))

    assert np.allclose(field.txx, 0.5)
    assert np.allclose(field.txy, 1.0)
    assert np.allclose(field.tyy, 1.5)


def test_case_01_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_01_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
