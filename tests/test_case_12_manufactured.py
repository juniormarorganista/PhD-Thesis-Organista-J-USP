import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_12_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_12_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_12_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters(beta=0.5, at=0.0))

    # Point (x, y) = (0, 0), with all stress coefficients equal to one.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 2.0)
    assert np.isclose(field.txy[0, 0], 2.0)
    assert np.isclose(field.tyy[0, 0], 2.0)


def test_case_12_known_center_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters(beta=0.5, at=0.0))

    center = (1, 1)  # y-index, x-index -> (x, y) = (0.5, 0.5)
    assert np.isclose(field.u[center], -0.25)
    assert np.isclose(field.v[center], 0.0, atol=1.0e-15)
    assert np.isclose(field.psi[center], -0.125)
    assert np.isclose(field.wz[center], 1.0 + np.pi**2 / 4.0)
    assert np.isclose(field.txx[center], 2.0)
    assert np.isclose(field.txy[center], 2.0)
    assert np.isclose(field.tyy[center], 2.0)


def test_case_12_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_12_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
