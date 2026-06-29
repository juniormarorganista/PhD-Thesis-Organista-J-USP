import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case06Parameters,
    evaluate_case_06_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_06_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_06_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_06_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_06_exact(grid, t=0.0, parameters=Case06Parameters(beta=0.5))

    # Point (x, y) = (0, 0): cos(0)=1, sin(0)=0, so u=0, v=0, wz=16*pi,
    # psi=-8/pi (decay=1). Stresses with all trigonometric coefficients
    # = 1 reduce to (1-beta)*(a1+a3)*(b1+b3) = 0.5*2*2 = 2.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 16.0 * np.pi)
    assert np.isclose(field.psi[0, 0], -8.0 / np.pi)
    assert np.isclose(field.txx[0, 0], 2.0)
    assert np.isclose(field.txy[0, 0], 2.0)
    assert np.isclose(field.tyy[0, 0], 2.0)


def test_case_06_velocity_fields_match_case_04() -> None:
    """Case 06's u/v/wz/psi are algebraically identical to Case 04's."""
    from mms_viscoelastic.manufactured_solutions.case_04 import evaluate_case_04_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_04 = evaluate_case_04_exact(grid, t=0.0)
    field_06 = evaluate_case_06_exact(grid, t=0.0)

    assert np.allclose(field_04.u, field_06.u)
    assert np.allclose(field_04.v, field_06.v)
    assert np.allclose(field_04.wz, field_06.wz)
    assert np.allclose(field_04.psi, field_06.psi)


def test_case_06_stress_fields_match_case_03() -> None:
    """Case 06's txx/txy/tyy are algebraically identical to Case 03's."""
    from mms_viscoelastic.manufactured_solutions.case_03 import evaluate_case_03_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_03 = evaluate_case_03_exact(grid, t=0.0)
    field_06 = evaluate_case_06_exact(grid, t=0.0)

    assert np.allclose(field_03.txx, field_06.txx)
    assert np.allclose(field_03.txy, field_06.txy)
    assert np.allclose(field_03.tyy, field_06.tyy)


def test_case_06_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_06_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
