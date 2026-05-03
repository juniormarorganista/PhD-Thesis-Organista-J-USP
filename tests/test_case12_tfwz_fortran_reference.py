"""Regression tests against Fortran-exported Case 12 `tfwz` reference data."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import compare_case12_tfwz_reference


REFERENCE_DIR = Path("validation") / "reference_data"


def test_case12_tfwz_matches_3x3_fortran_reference() -> None:
    comparison = compare_case12_tfwz_reference(
        REFERENCE_DIR / "case12_tfwz_3x3_default.csv"
    )

    assert comparison.passed
    assert comparison.nx == 3
    assert comparison.ny == 3
    assert comparison.max_abs_error < 1.0e-11


def test_case12_tfwz_matches_5x5_fortran_reference() -> None:
    comparison = compare_case12_tfwz_reference(
        REFERENCE_DIR / "case12_tfwz_5x5_default.csv"
    )

    assert comparison.passed
    assert comparison.nx == 5
    assert comparison.ny == 5
    assert comparison.max_abs_error < 1.0e-11
