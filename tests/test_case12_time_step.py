"""Regression tests for the Case 12 time-marching step against a Fortran reference CSV."""

from __future__ import annotations

from pathlib import Path

import pytest

from mms_viscoelastic.validation import (
    CASE12_TIME_STEP_COMPONENTS,
    compare_case12_time_step_reference,
)

REFERENCE_DIR = Path("validation") / "reference_data"


@pytest.mark.parametrize("component", CASE12_TIME_STEP_COMPONENTS)
def test_case12_time_step_reference_9x9(component: str) -> None:
    result = compare_case12_time_step_reference(
        REFERENCE_DIR / "case12_time_step_output_9x9_dt1e-05.csv",
        component=component,
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-12
    assert result.max_rel_error < 1.0e-12
