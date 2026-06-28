"""Regression tests against Fortran-exported Case 02 field/source-term reference data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.io.reference import read_reference_csv
from mms_viscoelastic.manufactured_solutions.case_02 import (
    Case02Parameters,
    evaluate_case_02_exact,
)
from mms_viscoelastic.source_terms.case_02 import (
    Case02SourceTermParameters,
    evaluate_case_02_source_terms,
)

REFERENCE_DIR = Path("validation") / "reference_data"


@pytest.mark.parametrize("size", [5, 9])
def test_case_02_fields_and_source_terms_match_fortran_reference(size: int) -> None:
    path = REFERENCE_DIR / f"case02_fields_{size}x{size}_default.csv"
    df = read_reference_csv(path)

    nx = ny = size
    assert nx * ny == len(df)

    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    field = evaluate_case_02_exact(grid, t=0.0, parameters=Case02Parameters())
    terms = evaluate_case_02_source_terms(grid, t=0.0, parameters=Case02SourceTermParameters())

    for column, predicted in (
        ("u", field.u),
        ("v", field.v),
        ("wz", field.wz),
        ("psi", field.psi),
        ("txx", field.txx),
        ("txy", field.txy),
        ("tyy", field.tyy),
        ("tfwz", terms.tfwz),
        ("tfTxx", terms.tfTxx),
        ("tfTxy", terms.tfTxy),
        ("tfTyy", terms.tfTyy),
    ):
        reference = df[column].to_numpy(dtype=float).reshape(ny, nx)
        max_abs_error = np.max(np.abs(predicted - reference))
        assert max_abs_error < 1.0e-10, f"{column}: max_abs_error={max_abs_error:.3e}"
