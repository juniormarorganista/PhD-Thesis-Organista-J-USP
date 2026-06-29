import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case11Parameters,
    evaluate_case_11_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_11_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_11_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_11_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_11_exact(grid, t=0.0, parameters=Case11Parameters(beta=0.5))

    # Point (x, y) = (0, 0): x_poly=(-1+x)*x vanishes at x=0, and sin(0)=0
    # also kills every velocity/vorticity/stream-function field. Stresses
    # with all polynomial coefficients = 1 reduce to (1-beta)*(1)*(1)=0.5.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 0.5)
    assert np.isclose(field.tyy[0, 0], 0.5)


def test_case_11_velocity_fields_match_case_10() -> None:
    """Case 11's u/v/wz/psi are algebraically identical to Case 10's."""
    from mms_viscoelastic.manufactured_solutions.case_10 import evaluate_case_10_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_10 = evaluate_case_10_exact(grid, t=0.0)
    field_11 = evaluate_case_11_exact(grid, t=0.0)

    assert np.allclose(field_10.u, field_11.u)
    assert np.allclose(field_10.v, field_11.v)
    assert np.allclose(field_10.wz, field_11.wz)
    assert np.allclose(field_10.psi, field_11.psi)


def test_case_11_stress_fields_match_case_02() -> None:
    """Case 11's txx/txy/tyy are algebraically identical to Case 02's."""
    from mms_viscoelastic.manufactured_solutions.case_02 import evaluate_case_02_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_02 = evaluate_case_02_exact(grid, t=0.0)
    field_11 = evaluate_case_11_exact(grid, t=0.0)

    assert np.allclose(field_02.txx, field_11.txx)
    assert np.allclose(field_02.txy, field_11.txy)
    assert np.allclose(field_02.tyy, field_11.tyy)


def test_case_11_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_11_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
