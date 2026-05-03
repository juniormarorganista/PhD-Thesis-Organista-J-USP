"""Contract tests for MMS Case 12 source terms."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import (
    Case12SourceTermParameters,
    available_case12_source_terms,
    case12_source_term_status,
    evaluate_case_12_source_terms,
    evaluate_case_12_tfTxx,
    evaluate_case_12_tfTxy,
    evaluate_case_12_tfwz,
    pending_case12_source_terms,
    ported_case12_source_terms,
)


def test_case12_source_term_registry() -> None:
    assert available_case12_source_terms() == ("tfwz", "tfTxx", "tfTxy", "tfTyy")
    assert ported_case12_source_terms() == ("tfwz", "tfTxx", "tfTxy")
    assert pending_case12_source_terms() == ("tfTyy",)


def test_case12_source_term_status_is_explicit() -> None:
    status = case12_source_term_status()

    assert set(status) == {"tfwz", "tfTxx", "tfTxy", "tfTyy"}
    assert "ported" in status["tfwz"]
    assert "ported" in status["tfTxx"]
    assert "ported" in status["tfTxy"]
    assert "pending" in status["tfTyy"]


def test_case12_source_term_parameters_defaults_are_safe() -> None:
    params = Case12SourceTermParameters()

    assert params.reynolds > 0.0
    assert params.wi > 0.0
    assert 0.0 <= params.beta <= 1.0


def test_case12_full_source_terms_still_raise_until_complete() -> None:
    with pytest.raises(NotImplementedError, match="not fully ported"):
        evaluate_case_12_source_terms()


def test_case12_ported_terms_shape_and_finite_values() -> None:
    grid = Grid2D.from_bounds(nx=9, ny=7)

    for evaluator in (
        evaluate_case_12_tfwz,
        evaluate_case_12_tfTxx,
        evaluate_case_12_tfTxy,
    ):
        values = evaluator(grid)
        assert values.shape == (7, 9)
        assert np.isfinite(values).all()


def test_case12_partial_source_term_object() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    terms = evaluate_case_12_source_terms(grid, allow_partial=True)

    assert terms.tfwz is not None
    assert terms.tfTxx is not None
    assert terms.tfTxy is not None
    assert terms.tfTyy is None

    assert terms.tfwz.shape == (4, 5)
    assert terms.tfTxx.shape == (4, 5)
    assert terms.tfTxy.shape == (4, 5)

    assert terms.available_components() == ("tfwz", "tfTxx", "tfTxy")
    assert terms.pending_components() == ("tfTyy",)


def test_case12_default_references_at_origin() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)

    assert np.isclose(evaluate_case_12_tfwz(grid)[0, 0], -4.0 * np.pi**2)
    assert np.isclose(evaluate_case_12_tfTxx(grid)[0, 0], 2.0)
    assert np.isclose(evaluate_case_12_tfTxy(grid)[0, 0], 2.0)
