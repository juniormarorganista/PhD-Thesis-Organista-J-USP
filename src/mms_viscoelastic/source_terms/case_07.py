"""Source terms for MMS Case 07.

This module ports the algebraic source terms from the original Fortran
routine:

``src/fortran_original/4Order_LPOG_Euler/case_07.f::source_term_case_07(t)``

Unlike every other ported case, case 07's velocity/vorticity fields decay
in time through a hyperbolic tangent (``tanh``, with ``1/cosh**2`` as its
time derivative) rather than ``exp(-at*t)``; the polymeric stress targets
``c1``, ``c2``, ``c3`` are the same spatially uniform constants used by
case 01.

This case also introduces ``dexp(8.0d0)``-style Fortran double-precision
literals (the ``d0`` exponent marker), which the mechanical
Fortran-to-Python script converts to the equivalent Python ``e0`` syntax.
Every term below was copied verbatim from that script's output.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE07_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case07SourceTermParameters:
    """Physical and constitutive parameters required by Case 07 source terms."""

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
class Case07SourceTerms:
    """Container for MMS Case 07 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def _case07_common(
    grid: Grid2D,
    t: float,
    parameters: Case07SourceTermParameters | None,
) -> dict[str, object]:
    """Collect common arrays and parameters for Case 07 source terms."""
    params = parameters or Case07SourceTermParameters()
    params.validate()

    return {
        "p": params,
        "x": grid.X,
        "y": grid.Y,
    }


def evaluate_case_07_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 07."""
    c = _case07_common(grid, t, parameters)
    p: Case07SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]

    Rey = p.reynolds
    betann = p.beta
    at = p.at

    return (
        -64
        * y
        * (
            2
            * at
            * Rey
            * (y**2 + 3 * (-1 + x) * x * ((-1 + x) * x + 2 * y**2))
            * (1.0e0 / np.cosh(4 - 8 * at * t)) ** 2
            + (-1 + np.tanh(4 - 8 * at * t))
            * (
                3 * betann * (1 + 6 * (-1 + x) * x + y**2)
                + 4
                * Rey
                * (-1 + x)
                * x
                * (-1 + 2 * x)
                * y**2
                * (y**2 + (-1 + x) * x * (-2 * (-1 + x) * x + 3 * y**2))
                + 4
                * Rey
                * (-1 + x)
                * x
                * (-1 + 2 * x)
                * y**2
                * (2 * (-1 + x) ** 2 * x**2 + (-1 - 3 * (-1 + x) * x) * y**2)
                * np.tanh(4 - 8 * at * t)
            )
        )
    ) / (3.0 * Rey)


def evaluate_case_07_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 07."""
    c = _case07_common(grid, t, parameters)
    p: Case07SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at

    return (
        (-1 + betann)
        * (
            -3 * alphaG * c2**2 * (np.exp(8.0e0) + np.exp(16 * at * t)) * Rey**2 * Wi
            - 3 * c1**2 * (np.exp(8.0e0) + np.exp(16 * at * t)) * (alphaG + epsylon) * Rey**2 * Wi
            - 3
            * c1
            * Rey
            * (
                c3 * np.exp(8.0e0) * epsylon * Rey * Wi
                + c3 * np.exp(16 * at * t) * epsylon * Rey * Wi
                + np.exp(8 + at * t) * (1 - at * Wi)
                + np.exp(17 * at * t) * (1 - at * Wi + 64 * Wi * (-1 + x) * x * (-1 + 2 * x) * (-1 + xi) * y**2)
            )
            + 32
            * np.exp(17 * at * t)
            * y
            * (
                6 * np.exp(at * t) * (-1 + x) * x * (-1 + 2 * x) * y
                + c2
                * Rey
                * Wi
                * (-3 * (-1 + x) ** 2 * x**2 * (-2 + xi) + (1 + 6 * (-1 + x) * x) * xi * y**2)
            )
        )
    ) / (3.0 * np.exp(2 * at * t) * (np.exp(8.0e0) + np.exp(16 * at * t)) * Rey)


def evaluate_case_07_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 07."""
    c = _case07_common(grid, t, parameters)
    p: Case07SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at

    return (
        (-1 + betann)
        * (
            3
            * c2
            * (np.exp(8.0e0) + np.exp(16 * at * t))
            * Rey
            * (-((c1 + c3) * (alphaG + epsylon) * Rey * Wi) + np.exp(at * t) * (-1 + at * Wi))
            + 16
            * np.exp(17 * at * t)
            * y
            * (
                2 * np.exp(at * t) * (3 * (-1 + x) ** 2 * x**2 + (-1 - 6 * (-1 + x) * x) * y**2)
                + Rey
                * Wi
                * (
                    -3 * (-1 + x) ** 2 * x**2 * (-2 * c3 + (c1 + c3) * xi)
                    + (1 + 6 * (-1 + x) * x) * (-2 * c1 + (c1 + c3) * xi) * y**2
                )
            )
        )
    ) / (3.0 * np.exp(2 * at * t) * (np.exp(8.0e0) + np.exp(16 * at * t)) * Rey)


def evaluate_case_07_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 07."""
    c = _case07_common(grid, t, parameters)
    p: Case07SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    c1, c2, c3 = p.c1, p.c2, p.c3
    at = p.at

    term_1 = (
        0.3333333333333333
        * (
            (-1 + betann)
            * (
                3 * alphaG * c2**2 * (np.exp(8.0e0) + np.exp(16 * at * t)) * Rey**2 * Wi
                + 3 * c3**2 * (np.exp(8.0e0) + np.exp(16 * at * t)) * (alphaG + epsylon) * Rey**2 * Wi
                + 32
                * np.exp(17 * at * t)
                * y
                * (
                    6 * np.exp(at * t) * (-1 + x) * x * (-1 + 2 * x) * y
                    + c2
                    * Rey
                    * Wi
                    * (3 * (-1 + x) ** 2 * x**2 * xi - (1 + 6 * (-1 + x) * x) * (-2 + xi) * y**2)
                )
                + 3
                * c3
                * Rey
                * (
                    c1 * np.exp(8.0e0) * epsylon * Rey * Wi
                    + c1 * np.exp(16 * at * t) * epsylon * Rey * Wi
                    + np.exp(8 + at * t) * (1 - at * Wi)
                    + np.exp(17 * at * t)
                    * (1 - at * Wi - 64 * Wi * (-1 + x) * x * (-1 + 2 * x) * (-1 + xi) * y**2)
                )
            )
        )
    ) / (np.exp(2 * at * t) * (np.exp(8.0e0) + np.exp(16 * at * t)) * Rey)

    return -term_1


def evaluate_case_07_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07SourceTermParameters | None = None,
) -> Case07SourceTerms:
    """Evaluate all Case 07 source terms."""
    return Case07SourceTerms(
        tfwz=evaluate_case_07_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_07_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_07_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_07_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE07_FORTRAN_ARRAYS",
    "Case07SourceTermParameters",
    "Case07SourceTerms",
    "evaluate_case_07_source_terms",
    "evaluate_case_07_tftxx",
    "evaluate_case_07_tftxy",
    "evaluate_case_07_tftyy",
    "evaluate_case_07_tfwz",
]
