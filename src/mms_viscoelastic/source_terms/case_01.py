"""Source terms for MMS Case 01.

This module ports the algebraic source terms from the original Fortran
routine:

``src/fortran_original/4Order_LPOG_Euler/case_01.f::source_term_case_01(t)``

Unlike case 12, case 01's manufactured stress field is spatially uniform
(constant target values ``c1``, ``c2``, ``c3`` from ``par.nn``), so the
source terms are correspondingly simpler algebraic expressions in ``x``,
``y`` (only ``tfwz`` depends on position; ``tfTxx``/``tfTxy``/``tfTyy`` mix
position-dependent terms from the velocity field with the constant stress
targets).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE01_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case01SourceTermParameters:
    """Physical and constitutive parameters required by Case 01 source terms."""

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
class Case01SourceTerms:
    """Container for MMS Case 01 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def _case01_common(
    grid: Grid2D,
    t: float,
    parameters: Case01SourceTermParameters | None,
) -> dict[str, object]:
    """Collect common arrays and parameters for Case 01 source terms."""
    params = parameters or Case01SourceTermParameters()
    params.validate()

    x = grid.X
    y = grid.Y

    return {
        "p": params,
        "x": x,
        "y": y,
        "exp_at": np.exp(params.at * t),
        "exp_2at": np.exp(2.0 * params.at * t),
    }


def evaluate_case_01_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 01."""
    c = _case01_common(grid, t, parameters)
    p: Case01SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    betann = p.beta
    at = p.at

    x_poly = (-1.0 + x) * x

    term_diffusion = 4.0 * betann * exp_at * (
        -1.0
        + 3.0 * y**2
        + 3.0
        * (
            (-2.0 + x) * (-1.0 + x) * x * (1.0 + x)
            + 12.0 * (-1.0 + x) * x * y**2
            + y**4
        )
    )

    term_convective = Rey * (
        at
        * exp_at
        * (
            -((-1.0 + x) ** 2 * x**2)
            + (-1.0 + 6.0 * (x - 2.0 * x**3 + x**4)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * y**4
        )
        + 32.0
        * (-1.0 + x)
        * x
        * (-1.0 + 2.0 * x)
        * y
        * (
            -((-1.0 + x) ** 2 * x**2)
            + (1.0 + x_poly) * (1.0 + 2.0 * x_poly) * y**2
            - 3.0 * (1.0 + x_poly) * (1.0 + 2.0 * x_poly) * y**4
            + 2.0 * (1.0 + 3.0 * x_poly) * y**6
        )
    )

    return 16.0 * (term_diffusion + term_convective) / (exp_2at * Rey)


def evaluate_case_01_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 01."""
    c = _case01_common(grid, t, parameters)
    p: Case01SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at
    c1 = p.c1
    c2 = p.c2
    c3 = p.c3

    velocity_term = (
        64.0 * exp_at * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y * (-1.0 + 2.0 * y**2) / Rey
    )

    relaxation_term = -c1 * (exp_at + (c1 + c3) * epsylon * Rey * Wi)

    giesekus_term = -(alphaG * (c1**2 + c2**2) * Rey * Wi)

    wi_inner = (
        at * c1 * exp_at
        + 64.0 * c1 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * (-1.0 + xi) * y * (1.0 - 2.0 * y**2)
        + 16.0
        * c2
        * (
            (-1.0 + x) ** 2 * x**2 * (-2.0 + xi)
            - (xi + 6.0 * (-1.0 + x) * x * ((-1.0 + x) * x * (-2.0 + xi) + xi)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * xi * y**4
        )
    )

    return (
        (-1.0 + betann)
        * (giesekus_term + relaxation_term + velocity_term + Wi * wi_inner)
        / exp_2at
    )


def evaluate_case_01_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 01."""
    c = _case01_common(grid, t, parameters)
    p: Case01SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at
    c1 = p.c1
    c2 = p.c2
    c3 = p.c3

    velocity_term = 8.0 * (
        2.0
        * exp_at
        * (
            -((-1.0 + x) ** 2 * x**2)
            + (1.0 + 6.0 * (-1.0 + x) * x * (1.0 + (-1.0 + x) * x)) * y**2
            + (-1.0 - 6.0 * (-1.0 + x) * x) * y**4
        )
        + Rey
        * Wi
        * (
            (-1.0 + x) ** 2 * x**2 * (c3 * (-2.0 + xi) + c1 * xi)
            - (
                2.0 * (-6.0 * c3 * (-1.0 + x) ** 2 * x**2 + c1 * (-1.0 - 6.0 * (-1.0 + x) * x))
                + (c1 + c3) * (1.0 + 6.0 * (-1.0 + x) * x * (1.0 + (-1.0 + x) * x)) * xi
            )
            * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * (c1 * (-2.0 + xi) + c3 * xi) * y**4
        )
    )

    relaxation_term = c2 * Rey * (
        -((c1 + c3) * (alphaG + epsylon) * Rey * Wi) + exp_at * (-1.0 + at * Wi)
    )

    return (-1.0 + betann) * (relaxation_term + velocity_term) / (exp_2at * Rey)


def evaluate_case_01_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 01."""
    c = _case01_common(grid, t, parameters)
    p: Case01SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at
    c1 = p.c1
    c2 = p.c2
    c3 = p.c3

    giesekus_term = -(alphaG * c2**2 * Rey**2 * Wi) - c3**2 * (alphaG + epsylon) * Rey**2 * Wi

    velocity_term = -64.0 * exp_at * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y * (-1.0 + 2.0 * y**2)

    c2_term = (
        16.0
        * c2
        * Rey
        * Wi
        * (
            (-1.0 + x) ** 2 * x**2 * xi
            - (-2.0 + xi + 6.0 * (-1.0 + x) * x * (-2.0 + xi + (-1.0 + x) * x * xi)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * (-2.0 + xi) * y**4
        )
    )

    c3_term = c3 * Rey * (
        -(c1 * epsylon * Rey * Wi)
        + exp_at * (-1.0 + at * Wi)
        + 64.0 * Wi * (-1.0 + x) * x * (-1.0 + 2.0 * x) * (-1.0 + xi) * y * (-1.0 + 2.0 * y**2)
    )

    return (
        (-1.0 + betann)
        * (giesekus_term + velocity_term + c2_term + c3_term)
        / (exp_2at * Rey)
    )


def evaluate_case_01_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01SourceTermParameters | None = None,
) -> Case01SourceTerms:
    """Evaluate all Case 01 source terms."""
    return Case01SourceTerms(
        tfwz=evaluate_case_01_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_01_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_01_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_01_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE01_FORTRAN_ARRAYS",
    "Case01SourceTermParameters",
    "Case01SourceTerms",
    "evaluate_case_01_source_terms",
    "evaluate_case_01_tftxx",
    "evaluate_case_01_tftxy",
    "evaluate_case_01_tftyy",
    "evaluate_case_01_tfwz",
]
