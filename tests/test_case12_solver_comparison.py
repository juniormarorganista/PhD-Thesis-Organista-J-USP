"""Regression tests for Case 12 t=0 exact fields against Fortran reference CSVs."""

from __future__ import annotations

from pathlib import Path

import pytest

from mms_viscoelastic.validation import CASE12_FIELD_COMPONENTS, compare_case12_field_reference

REFERENCE_DIR = Path("validation") / "reference_data"


@pytest.mark.parametrize("component", CASE12_FIELD_COMPONENTS)
def test_case12_field_reference_3x3(component: str) -> None:
    result = compare_case12_field_reference(
        REFERENCE_DIR / "case12_fields_3x3_default.csv",
        component=component,
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-10
    assert result.max_rel_error < 1.0e-10


@pytest.mark.parametrize("component", CASE12_FIELD_COMPONENTS)
def test_case12_field_reference_5x5(component: str) -> None:
    result = compare_case12_field_reference(
        REFERENCE_DIR / "case12_fields_5x5_default.csv",
        component=component,
    )

    assert result.passed
    assert result.max_abs_error < 1.0e-10
    assert result.max_rel_error < 1.0e-10
