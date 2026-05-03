"""Regression tests for reference-comparison grid metadata."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import (
    compare_case12_tftxx_reference,
    compare_case12_tftxy_reference,
    compare_case12_tfwz_reference,
)


REFERENCE_DIR = Path("validation") / "reference_data"


def test_tfwz_comparison_exposes_grid_shape() -> None:
    comparison = compare_case12_tfwz_reference(
        REFERENCE_DIR / "case12_tfwz_3x3_default.csv"
    )

    assert comparison.nx == 3
    assert comparison.ny == 3
    assert comparison.n_points == 9


def test_tftxx_comparison_exposes_grid_shape() -> None:
    comparison = compare_case12_tftxx_reference(
        REFERENCE_DIR / "case12_tfTxx_5x5_default.csv"
    )

    assert comparison.nx == 5
    assert comparison.ny == 5
    assert comparison.n_points == 25


def test_tftxy_comparison_exposes_grid_shape() -> None:
    comparison = compare_case12_tftxy_reference(
        REFERENCE_DIR / "case12_tfTxy_5x5_default.csv"
    )

    assert comparison.nx == 5
    assert comparison.ny == 5
    assert comparison.n_points == 25
