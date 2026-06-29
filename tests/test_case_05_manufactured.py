import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case05Parameters,
    evaluate_case_05_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_05_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_05_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_05_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_05_exact(grid, t=0.0, parameters=Case05Parameters(beta=0.5))

    # Point (x, y) = (0, 0): cos(0)=1, sin(0)=0, so u=0, v=0, wz=16*pi,
    # psi=-8/pi (decay=1). Stresses with all polynomial coefficients = 1
    # reduce to (1-beta)*(1)*(1) = 0.5.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 16.0 * np.pi)
    assert np.isclose(field.psi[0, 0], -8.0 / np.pi)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 0.5)
    assert np.isclose(field.tyy[0, 0], 0.5)


def test_case_05_velocity_fields_match_case_04() -> None:
    """Case 05's u/v/wz/psi are algebraically identical to Case 04's."""
    from mms_viscoelastic.manufactured_solutions.case_04 import evaluate_case_04_exact

    grid = Grid2D.from_bounds(nx=9, ny=9)
    field_04 = evaluate_case_04_exact(grid, t=0.0)
    field_05 = evaluate_case_05_exact(grid, t=0.0)

    assert np.allclose(field_04.u, field_05.u)
    assert np.allclose(field_04.v, field_05.v)
    assert np.allclose(field_04.wz, field_05.wz)
    assert np.allclose(field_04.psi, field_05.psi)


def test_case_05_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_05_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
