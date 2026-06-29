import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case03Parameters,
    evaluate_case_03_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_03_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_03_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_03_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_03_exact(grid, t=0.0, parameters=Case03Parameters(beta=0.5))

    # Point (x, y) = (0, 0): velocity/vorticity/stream-function vanish
    # (factor x or y^2*(-1+y^2) is zero). Stresses with all trigonometric
    # coefficients = 1 reduce to (1-beta)*(a1+a3)*(b1+b3) = 0.5*2*2 = 2.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 2.0)
    assert np.isclose(field.txy[0, 0], 2.0)
    assert np.isclose(field.tyy[0, 0], 2.0)


def test_case_03_velocity_fields_match_case_01() -> None:
    """Case 03's u/v/wz/psi are algebraically identical to Case 01's."""
    from mms_viscoelastic.manufactured_solutions.case_01 import evaluate_case_01_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_01 = evaluate_case_01_exact(grid, t=0.0)
    field_03 = evaluate_case_03_exact(grid, t=0.0)

    assert np.allclose(field_01.u, field_03.u)
    assert np.allclose(field_01.v, field_03.v)
    assert np.allclose(field_01.wz, field_03.wz)
    assert np.allclose(field_01.psi, field_03.psi)


def test_case_03_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_03_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
