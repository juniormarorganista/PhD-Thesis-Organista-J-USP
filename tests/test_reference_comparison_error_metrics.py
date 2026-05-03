"""Regression tests for reference-comparison error metrics."""

from __future__ import annotations

from pathlib import Path

from mms_viscoelastic.validation import compare_case12_tfwz_reference


REFERENCE_DIR = Path("validation") / "reference_data"


def test_reference_comparison_exposes_mean_errors() -> None:
    comparison = compare_case12_tfwz_reference(
        REFERENCE_DIR / "case12_tfwz_3x3_default.csv"
    )

    assert comparison.nx == 3
    assert comparison.ny == 3
    assert comparison.n_points == 9
    assert comparison.max_abs_error >= 0.0
    assert comparison.mean_abs_error >= 0.0
    assert comparison.max_rel_error >= 0.0
    assert comparison.mean_rel_error >= 0.0
    assert comparison.mean_abs_error <= comparison.max_abs_error
    assert comparison.mean_rel_error <= comparison.max_rel_error
