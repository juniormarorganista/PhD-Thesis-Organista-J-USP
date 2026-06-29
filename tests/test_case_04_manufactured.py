import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case04Parameters,
    evaluate_case_04_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_04_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_04_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_04_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_04_exact(grid, t=0.0, parameters=Case04Parameters(beta=0.5))

    # Point (x, y) = (0, 0): cos(0)=1, sin(0)=0, so u=0, v=0, wz=16*pi,
    # psi=-8/pi (decay=1). Stresses with defaults c1=1, c2=2, c3=3 reduce
    # to (1-beta)*c_i = 0.5*c_i.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 16.0 * np.pi)
    assert np.isclose(field.psi[0, 0], -8.0 / np.pi)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 1.0)
    assert np.isclose(field.tyy[0, 0], 1.5)


def test_case_04_stress_fields_match_case_01() -> None:
    """Case 04's txx/txy/tyy are algebraically identical to Case 01's."""
    from mms_viscoelastic.manufactured_solutions.case_01 import evaluate_case_01_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_01 = evaluate_case_01_exact(grid, t=0.0)
    field_04 = evaluate_case_04_exact(grid, t=0.0)

    assert np.allclose(field_01.txx, field_04.txx)
    assert np.allclose(field_01.txy, field_04.txy)
    assert np.allclose(field_01.tyy, field_04.tyy)


def test_case_04_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_04_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
