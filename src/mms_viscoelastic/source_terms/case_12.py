"""Source-term contract for MMS Case 12.

This module defines the public Python contract for the manufactured source
terms of Case 12. The full algebraic port from the original Fortran routine
``source_term_case_12(t)`` is intentionally not faked here.

The original Fortran routine fills four arrays:

- ``tfwz``: source term for the vorticity equation;
- ``tfTxx``: source term for the ``Txx`` constitutive equation;
- ``tfTxy``: source term for the ``Txy`` constitutive equation;
- ``tfTyy``: source term for the ``Tyy`` constitutive equation.

Until the exact algebraic expressions are fully ported and validated, calling
``evaluate_case_12_source_terms`` raises ``NotImplementedError``. This is safer
than returning zeros or approximate placeholders that could contaminate later
verification or machine-learning datasets.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


CASE12_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case12SourceTermParameters:
    """Physical and constitutive parameters required by Case 12 source terms."""

    reynolds: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    at: float = 0.0

    # Coefficients used in the manufactured tensor fields.
    a1xx: float = 1.0
    a2xx: float = 1.0
    a3xx: float = 1.0
    b1xx: float = 1.0
    b2xx: float = 1.0
    b3xx: float = 1.0

    a1xy: float = 1.0
    a2xy: float = 1.0
    a3xy: float = 1.0
    b1xy: float = 1.0
    b2xy: float = 1.0
    b3xy: float = 1.0

    a1yy: float = 1.0
    a2yy: float = 1.0
    a3yy: float = 1.0
    b1yy: float = 1.0
    b2yy: float = 1.0
    b3yy: float = 1.0


@dataclass(frozen=True)
class Case12SourceTerms:
    """Container for MMS Case 12 source-term arrays."""

    tfwz: np.ndarray
    tfTxx: np.ndarray
    tfTxy: np.ndarray
    tfTyy: np.ndarray


def available_case12_source_terms() -> tuple[str, ...]:
    """Return the Fortran source-term arrays targeted by the Python port."""

    return CASE12_FORTRAN_ARRAYS


def case12_source_term_status() -> dict[str, str]:
    """Return current implementation status for each source-term component."""

    return {
        "tfwz": "mapped from Fortran, algebraic Python port pending",
        "tfTxx": "mapped from Fortran, algebraic Python port pending",
        "tfTxy": "mapped from Fortran, algebraic Python port pending",
        "tfTyy": "mapped from Fortran, algebraic Python port pending",
    }


def evaluate_case_12_source_terms(
    *args: object,
    **kwargs: object,
) -> Case12SourceTerms:
    """Evaluate MMS Case 12 source terms.

    The function currently defines the public contract only. The exact
    expressions from ``source_term_case_12(t)`` still need to be ported and
    independently validated against the original Fortran code.

    Raises
    ------
    NotImplementedError
        Always, until the full source-term expressions are ported.
    """

    raise NotImplementedError(
        "MMS Case 12 source terms are mapped but not yet ported. "
        "Port the exact formulas from "
        "src/fortran_original/4Order_LPOG_Euler/case_12.f::source_term_case_12 "
        "before using this function in verification or ML workflows."
    )


__all__ = [
    "CASE12_FORTRAN_ARRAYS",
    "Case12SourceTermParameters",
    "Case12SourceTerms",
    "available_case12_source_terms",
    "case12_source_term_status",
    "evaluate_case_12_source_terms",
]
