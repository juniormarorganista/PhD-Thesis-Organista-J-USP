import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case10Parameters,
    evaluate_case_10_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_10_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_10_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_10_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_10_exact(grid, t=0.0, parameters=Case10Parameters(beta=0.5))

    # Point (x, y) = (0, 0): x_poly=(-1+x)*x vanishes at x=0, and sin(0)=0
    # also kills every velocity/vorticity/stream-function field. Stresses
    # with defaults c1=1, c2=2, c3=3 reduce to (1-beta)*c_i = 0.5*c_i.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 1.0)
    assert np.isclose(field.tyy[0, 0], 1.5)


def test_case_10_stress_fields_match_case_01() -> None:
    """Case 10's txx/txy/tyy are algebraically identical to Case 01's."""
    from mms_viscoelastic.manufactured_solutions.case_01 import evaluate_case_01_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_01 = evaluate_case_01_exact(grid, t=0.0)
    field_10 = evaluate_case_10_exact(grid, t=0.0)

    assert np.allclose(field_01.txx, field_10.txx)
    assert np.allclose(field_01.txy, field_10.txy)
    assert np.allclose(field_01.tyy, field_10.tyy)


def test_case_10_velocity_fields_match_case_12() -> None:
    """Case 10's u/v/wz/psi are algebraically identical to Case 12's."""
    from mms_viscoelastic.manufactured_solutions.case_12 import evaluate_case_12_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_12 = evaluate_case_12_exact(grid, t=0.0)
    field_10 = evaluate_case_10_exact(grid, t=0.0)

    assert np.allclose(field_12.u, field_10.u)
    assert np.allclose(field_12.v, field_10.v)
    assert np.allclose(field_12.wz, field_10.wz)
    assert np.allclose(field_12.psi, field_10.psi)


def test_case_10_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_10_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
