"""Regression tests for validation-grid reconstruction."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import (
    compare_case12_tftxx_reference,
    compare_case12_tftxy_reference,
    compare_case12_tfwz_reference,
)


REFERENCE_DIR = Path("validation") / "reference_data"


def test_validation_grid_api_tfwz_reference_runs() -> None:
    result = compare_case12_tfwz_reference(
        REFERENCE_DIR / "case12_tfwz_3x3_default.csv"
    )

    assert result.passed


def test_validation_grid_api_tftxx_reference_runs() -> None:
    result = compare_case12_tftxx_reference(
        REFERENCE_DIR / "case12_tfTxx_3x3_default.csv"
    )

    assert result.passed


def test_validation_grid_api_tftxy_reference_runs() -> None:
    result = compare_case12_tftxy_reference(
        REFERENCE_DIR / "case12_tfTxy_3x3_default.csv"
    )

    assert result.passed
