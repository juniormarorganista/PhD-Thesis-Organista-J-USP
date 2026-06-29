import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_14_exact
from mms_viscoelastic.ml import field_to_feature_table
from mms_viscoelastic.source_terms import evaluate_case_14_source_terms


def test_case_14_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_14_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_14_fields_are_identically_zero() -> None:
    """Case 14 is not an MMS case: every exact field is identically zero
    (it is a smooth-startup lid-driven cavity flow benchmark)."""
    grid = Grid2D.from_bounds(nx=9, ny=9)
    field = evaluate_case_14_exact(grid, t=1.5)

    assert np.all(field.u == 0.0)
    assert np.all(field.v == 0.0)
    assert np.all(field.wz == 0.0)
    assert np.all(field.psi == 0.0)
    assert np.all(field.txx == 0.0)
    assert np.all(field.txy == 0.0)
    assert np.all(field.tyy == 0.0)


def test_case_14_source_terms_are_identically_zero() -> None:
    grid = Grid2D.from_bounds(nx=9, ny=9)
    terms = evaluate_case_14_source_terms(grid, t=1.5)

    assert np.all(terms.tfwz == 0.0)
    assert np.all(terms.tfTxx == 0.0)
    assert np.all(terms.tfTxy == 0.0)
    assert np.all(terms.tfTyy == 0.0)


def test_case_14_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_14_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
