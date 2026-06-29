"""Source terms for MMS Case 14.

Unlike cases 01-12, case 14 is **not** a method-of-manufactured-solutions
verification case: it is a smooth-startup variant of the lid-driven
cavity flow benchmark. Its ``source_term_case_14(t)`` subroutine sets
every source term (``tfwz``, ``tfTxx``, ``tfTxy``, ``tfTyy``) to ``0.0``
identically -- there is no forcing to derive, which is also why no
Mathematica notebook exists for this case.

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

CASE14_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case14SourceTermParameters:
    """Parameters for case 14's source terms. Unused: every term is identically zero."""

    def validate(self) -> None:
        """No parameters to validate; present for API parity with other cases."""
        return


@dataclass(frozen=True)
class Case14SourceTerms:
    """Container for MMS Case 14 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def evaluate_case_14_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 14 (identically zero)."""
    (parameters or Case14SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_14_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 14 (identically zero)."""
    (parameters or Case14SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_14_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 14 (identically zero)."""
    (parameters or Case14SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_14_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 14 (identically zero)."""
    (parameters or Case14SourceTermParameters()).validate()
    return np.zeros((grid.ny, grid.nx))


def evaluate_case_14_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14SourceTermParameters | None = None,
) -> Case14SourceTerms:
    """Evaluate all Case 14 source terms."""
    return Case14SourceTerms(
        tfwz=evaluate_case_14_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_14_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_14_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_14_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE14_FORTRAN_ARRAYS",
    "Case14SourceTermParameters",
    "Case14SourceTerms",
    "evaluate_case_14_source_terms",
    "evaluate_case_14_tftxx",
    "evaluate_case_14_tftxy",
    "evaluate_case_14_tftyy",
    "evaluate_case_14_tfwz",
]
