"""Source terms for MMS Case 13.

Unlike cases 01-12, case 13 is **not** a method-of-manufactured-solutions
verification case: it is the classic lid-driven cavity flow benchmark.
Its ``source_term_case_13(t)`` subroutine sets every source term
(``tfwz``, ``tfTxx``, ``tfTxy``, ``tfTyy``) to ``0.0`` identically --
there is no forcing to derive, which is also why no Mathematica notebook
exists for this case.

This module is ported only for structural parity with the other 12
cases; it has no verification value beyond confirming the source terms
are identically zero.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE13_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case13SourceTermParameters:
    """Parameters for case 13's source terms. Unused: every term is identically zero."""

    def validate(self) -> None:
        """No parameters to validate; present for API parity with other cases."""
        return


@dataclass(frozen=True)
class Case13SourceTerms:
    """Container for MMS Case 13 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def evaluate_case_13_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 13 (identically zero)."""
    (parameters or Case13SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_13_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 13 (identically zero)."""
    (parameters or Case13SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_13_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 13 (identically zero)."""
    (parameters or Case13SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_13_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 13 (identically zero)."""
    (parameters or Case13SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_13_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13SourceTermParameters | None = None,
) -> Case13SourceTerms:
    """Evaluate all Case 13 source terms."""
    return Case13SourceTerms(
        tfwz=evaluate_case_13_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_13_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_13_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_13_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE13_FORTRAN_ARRAYS",
    "Case13SourceTermParameters",
    "Case13SourceTerms",
    "evaluate_case_13_source_terms",
    "evaluate_case_13_tftxx",
    "evaluate_case_13_tftxy",
    "evaluate_case_13_tftyy",
    "evaluate_case_13_tfwz",
]
