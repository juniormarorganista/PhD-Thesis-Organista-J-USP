"""Contract tests for MMS Case 12 source terms."""

from __future__ import annotations

import pytest

from mms_viscoelastic.source_terms import (
    Case12SourceTermParameters,
    available_case12_source_terms,
    case12_source_term_status,
    evaluate_case_12_source_terms,
)


def test_case12_source_term_registry() -> None:
    assert available_case12_source_terms() == ("tfwz", "tfTxx", "tfTxy", "tfTyy")


def test_case12_source_term_status_is_explicit() -> None:
    status = case12_source_term_status()

    assert set(status) == {"tfwz", "tfTxx", "tfTxy", "tfTyy"}
    assert all("pending" in value for value in status.values())


def test_case12_source_term_parameters_defaults_are_safe() -> None:
    params = Case12SourceTermParameters()

    assert params.reynolds > 0.0
    assert params.wi > 0.0
    assert 0.0 <= params.beta <= 1.0


def test_case12_source_terms_raise_until_ported() -> None:
    with pytest.raises(NotImplementedError, match="not yet ported"):
        evaluate_case_12_source_terms()
