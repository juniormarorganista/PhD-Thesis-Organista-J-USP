import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case08Parameters,
    evaluate_case_08_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_08_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_08_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_08_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_08_exact(grid, t=0.0, parameters=Case08Parameters(beta=0.5))

    # Point (x, y) = (0, 0): the x_poly=(1-x)^2*x^2 factor vanishes at x=0
    # for u/wz/psi, and y=0 kills v/wz/psi too -- all velocity/vorticity/
    # stream-function fields are exactly zero there for any t. Stresses
    # with all polynomial coefficients = 1 reduce to (1-beta)*(1)*(1)=0.5.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 0.5)
    assert np.isclose(field.txy[0, 0], 0.5)
    assert np.isclose(field.tyy[0, 0], 0.5)


def test_case_08_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_08_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
