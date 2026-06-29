"""Source terms for MMS Case 04.

This module ports the algebraic source terms from the original Fortran
routine:

``src/fortran_original/4Order_LPOG_Euler/case_04.f::source_term_case_04(t)``

The velocity/vorticity fields follow a Taylor-Green-vortex trigonometric
ansatz (distinct from cases 01/02/03 and case 12); the polymeric stress
targets ``c1``, ``c2``, ``c3`` are the same spatially uniform constants
used by case 01.

These source terms are considerably more compact than cases 02/03 (the
largest, ``tfTxx``/``tfTyy``, span 19 matched parenthesis pairs). Because
they read as short, the outer ``(-1 + betann)`` factor multiplying
``tfTxx``/``tfTxy``/``tfTyy`` was retyped by hand as ``(1 - betann)``
instead of copied verbatim from the mechanically-extracted expression --
a sign flip, since ``-1 + betann`` and ``1 - betann`` are negatives of
each other. The Fortran cross-reference caught it immediately (a
constant ~100-200 absolute error at every grid point, independent of
``x``/``y``, which is the signature of a sign error in a term that does
not depend on position). Fixed by restoring ``(-1.0 + betann)`` to match
the original. The lesson generalizes: "short enough to read at a glance"
is not a reason to skip the verbatim-copy discipline used for cases 02/03.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE04_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case04SourceTermParameters:
    """Physical and constitutive parameters required by Case 04 source terms."""

    reynolds: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    at: float = 0.0
    c1: float = 1.0
    c2: float = 2.0
    c3: float = 3.0

    def validate(self) -> None:
        """Validate parameters needed by the ported source terms."""
        if self.reynolds == 0.0:
            raise ValueError("reynolds must be nonzero.")
        if self.wi == 0.0:
            raise ValueError("wi must be nonzero.")
        if not np.isfinite(self.beta):
            raise ValueError("beta must be finite.")


@dataclass(frozen=True)
class Case04SourceTerms:
    """Container for MMS Case 04 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def _case04_common(
    grid: Grid2D,
    t: float,
    parameters: Case04SourceTermParameters | None,
) -> dict[str, object]:
    """Collect common arrays and parameters for Case 04 source terms."""
    params = parameters or Case04SourceTermParameters()
    params.validate()

    x = grid.X
    y = grid.Y
    pi = np.pi

    return {
        "p": params,
        "x": x,
        "y": y,
        "cos_px": np.cos(pi * x),
        "sin_px": np.sin(pi * x),
        "cos_py": np.cos(pi * y),
        "sin_py": np.sin(pi * y),
        "exp_at": np.exp(params.at * t),
        "exp_2at": np.exp(2.0 * params.at * t),
    }


def evaluate_case_04_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 04."""
    c = _case04_common(grid, t, parameters)
    p: Case04SourceTermParameters = c["p"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]

    Rey = p.reynolds
    betann = p.beta
    at = p.at
    pi = np.pi

    return (16.0 * pi * (-2.0 * betann * pi**2 + at * Rey) * cos_px * cos_py) / (exp_at * Rey)


def evaluate_case_04_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 04."""
    c = _case04_common(grid, t, parameters)
    p: Case04SourceTermParameters = c["p"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at
    pi = np.pi

    return (
        (-1.0 + betann)
        * (
            Rey * (-((alphaG * (c1**2 + c2**2) + c1 * (c1 + c3) * epsylon) * Rey * Wi) + c1 * exp_at * (-1.0 + at * Wi))
            + 16.0 * c2 * pi * Rey * Wi * cos_px * cos_py
            - 16.0 * pi * (exp_at - c1 * Rey * Wi * (-1.0 + xi)) * sin_px * sin_py
        )
    ) / (exp_2at * Rey)


def evaluate_case_04_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 04."""
    c = _case04_common(grid, t, parameters)
    p: Case04SourceTermParameters = c["p"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at
    pi = np.pi

    return (
        (-1.0 + betann)
        * (
            -(c2 * (c1 + c3) * (alphaG + epsylon) * Rey * Wi)
            + c2 * exp_at * (-1.0 + at * Wi)
            + 8.0 * (-c1 + c3) * pi * Wi * cos_px * cos_py
        )
    ) / exp_2at


def evaluate_case_04_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 04."""
    c = _case04_common(grid, t, parameters)
    p: Case04SourceTermParameters = c["p"]  # type: ignore[assignment]
    cos_px: FloatArray = c["cos_px"]  # type: ignore[assignment]
    sin_px: FloatArray = c["sin_px"]  # type: ignore[assignment]
    cos_py: FloatArray = c["cos_py"]  # type: ignore[assignment]
    sin_py: FloatArray = c["sin_py"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at
    pi = np.pi

    return (
        (-1.0 + betann)
        * (
            Rey * (-((alphaG * (c2**2 + c3**2) + c3 * (c1 + c3) * epsylon) * Rey * Wi) + c3 * exp_at * (-1.0 + at * Wi))
            - 16.0 * c2 * pi * Rey * Wi * cos_px * cos_py
            + 16.0 * pi * (exp_at - c3 * Rey * Wi * (-1.0 + xi)) * sin_px * sin_py
        )
    ) / (exp_2at * Rey)


def evaluate_case_04_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04SourceTermParameters | None = None,
) -> Case04SourceTerms:
    """Evaluate all Case 04 source terms."""
    return Case04SourceTerms(
        tfwz=evaluate_case_04_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_04_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_04_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_04_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE04_FORTRAN_ARRAYS",
    "Case04SourceTermParameters",
    "Case04SourceTerms",
    "evaluate_case_04_source_terms",
    "evaluate_case_04_tftxx",
    "evaluate_case_04_tftxy",
    "evaluate_case_04_tftyy",
    "evaluate_case_04_tfwz",
]
