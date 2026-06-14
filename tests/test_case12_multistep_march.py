"""Regression tests for multi-step Case 12 time marching against a Fortran reference CSV."""

from __future__ import annotations

from pathlib import Path

import pytest

from mms_viscoelastic.validation import (
    CASE12_TIME_STEP_COMPONENTS,
    compare_case12_multistep_reference,
)

REFERENCE_DIR = Path("validation") / "reference_data"
N_STEPS = 10


@pytest.mark.parametrize("component", CASE12_TIME_STEP_COMPONENTS)
def test_case12_multistep_reference_9x9(component: str) -> None:
    result = compare_case12_multistep_reference(
        REFERENCE_DIR / "case12_multistep_output_9x9_dt1e-05_n10.csv",
        n_steps=N_STEPS,
        component=component,
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-12
    assert result.max_rel_error < 1.0e-12
