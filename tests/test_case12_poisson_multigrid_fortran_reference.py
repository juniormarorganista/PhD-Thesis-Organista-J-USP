"""Regression test for the 4th-order multigrid Poisson solver against a Fortran reference CSV."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import compare_case12_poisson_multigrid_reference

REFERENCE_DIR = Path("validation") / "reference_data"


def test_case12_poisson_multigrid_reference_33x33() -> None:
    result = compare_case12_poisson_multigrid_reference(
        REFERENCE_DIR / "case12_poisson_multigrid_output_33x33.csv",
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-12
    assert result.max_rel_error < 1.0e-12
