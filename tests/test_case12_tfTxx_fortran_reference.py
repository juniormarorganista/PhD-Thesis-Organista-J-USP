"""Regression tests for Case 12 ``tfTxx`` against Fortran reference CSV files."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import compare_case12_tftxx_reference


REFERENCE_DIR = Path("validation") / "reference_data"


def test_case12_tftxx_reference_3x3() -> None:
    result = compare_case12_tftxx_reference(
        REFERENCE_DIR / "case12_tfTxx_3x3_default.csv"
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-10
    assert result.max_rel_error < 1.0e-10


def test_case12_tftxx_reference_5x5() -> None:
    result = compare_case12_tftxx_reference(
        REFERENCE_DIR / "case12_tfTxx_5x5_default.csv"
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-10
    assert result.max_rel_error < 1.0e-10
