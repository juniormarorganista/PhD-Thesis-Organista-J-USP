"""Source terms for MMS Case 02.

This module ports the algebraic source terms from the original Fortran
routine:

``src/fortran_original/4Order_LPOG_Euler/case_02.f::source_term_case_02(t)``

The polymeric tensor shape factors ``Axx = a1xx + a2xx x + a3xx x^2``,
``Bxx = b1xx + b2xx y + b3xx y^2`` (and the ``xy``/``yy`` analogues) follow
the same naming as case 12's trigonometric shape factors, but are
quadratic polynomials here instead of trigonometric functions.

Relaxation-term simplification
-------------------------------
Each tensor source term contains a sub-expression of the form::

    a1xx*Bxx + a1yy*Byy
    + x*(a2yy*b1yy + a2yy*y*(b2yy + b3yy*y)
         + a2xx*Bxx + a3xx*x*Bxx + a3yy*x*Byy)

which is algebraically identical to ``Axx*Bxx + Ayy*Byy`` (verified with
SymPy; the two differ by exactly ``0`` after expansion). The Python port
uses this simplified, already-factored form -- the same one used for the
relaxation term in case 12's source terms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE02_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case02SourceTermParameters:
    """Physical and constitutive parameters required by Case 02 source terms."""

    reynolds: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    at: float = 0.0

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

    def validate(self) -> None:
        """Validate parameters needed by the ported source terms."""
        if self.reynolds == 0.0:
            raise ValueError("reynolds must be nonzero.")
        if self.wi == 0.0:
            raise ValueError("wi must be nonzero.")
        if not np.isfinite(self.beta):
            raise ValueError("beta must be finite.")


@dataclass(frozen=True)
class Case02SourceTerms:
    """Container for MMS Case 02 source-term arrays."""

    tfwz: FloatArray
    tfTxx: FloatArray
    tfTxy: FloatArray
    tfTyy: FloatArray


def _case02_common(
    grid: Grid2D,
    t: float,
    parameters: Case02SourceTermParameters | None,
) -> dict[str, object]:
    """Collect common arrays and parameters for Case 02 source terms."""
    params = parameters or Case02SourceTermParameters()
    params.validate()

    x = grid.X
    y = grid.Y

    Axx = params.a1xx + x * (params.a2xx + params.a3xx * x)
    Bxx = params.b1xx + y * (params.b2xx + params.b3xx * y)
    Axy = params.a1xy + x * (params.a2xy + params.a3xy * x)
    Bxy = params.b1xy + y * (params.b2xy + params.b3xy * y)
    Ayy = params.a1yy + x * (params.a2yy + params.a3yy * x)
    Byy = params.b1yy + y * (params.b2yy + params.b3yy * y)

    return {
        "p": params,
        "x": x,
        "y": y,
        "exp_at": np.exp(params.at * t),
        "exp_2at": np.exp(2.0 * params.at * t),
        "Axx": Axx,
        "Bxx": Bxx,
        "Axy": Axy,
        "Bxy": Bxy,
        "Ayy": Ayy,
        "Byy": Byy,
        "tensor_sum": Axx * Bxx + Ayy * Byy,
    }


def evaluate_case_02_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 02."""
    c = _case02_common(grid, t, parameters)
    p: Case02SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]

    Rey = p.reynolds
    betann = p.beta
    at = p.at
    a2xx, a3xx, b2xx, b3xx = p.a2xx, p.a3xx, p.b2xx, p.b3xx
    a1xy, a2xy, a3xy, b1xy, b2xy, b3xy = p.a1xy, p.a2xy, p.a3xy, p.b1xy, p.b2xy, p.b3xy
    a2yy, a3yy, b2yy, b3yy = p.a2yy, p.a3yy, p.b2yy, p.b3yy

    betann_term = (
        -64.0
        + a2yy * b2yy * Rey
        - 2.0 * a1xy * b3xy * Rey
        + 384.0 * x
        - 2.0 * a3xx * b2xx * Rey * x
        + 2.0 * a3yy * b2yy * Rey * x
        - 2.0 * a2xy * b3xy * Rey * x
        - 192.0 * x**2
        - 384.0 * x**3
        + 192.0 * x**4
        + 2.0 * Rey * (a2yy * b3yy - 2.0 * a3xx * b3xx * x + 2.0 * a3yy * b3yy * x) * y
        + 192.0 * (1.0 + 12.0 * (-1.0 + x) * x) * y**2
        + 192.0 * y**4
        - a2xx * Rey * (b2xx + 2.0 * b3xx * y)
        + 2.0 * a3xy * Rey * (b1xy - b3xy * x**2 + b2xy * y + b3xy * y**2)
    )

    convective_term = (
        -512.0
        * (-1.0 + x)
        * x
        * (-1.0 + 2.0 * x)
        * y
        * (
            (-1.0 + x) ** 2 * x**2
            - (1.0 + (-1.0 + x) * x) * (1.0 + 2.0 * (-1.0 + x) * x) * y**2
            + 3.0 * (1.0 + (-1.0 + x) * x) * (1.0 + 2.0 * (-1.0 + x) * x) * y**4
            - 2.0 * (1.0 + 3.0 * (-1.0 + x) * x) * y**6
        )
    )

    exp_term = (
        -(a2yy * b2yy)
        + 2.0 * a1xy * b3xy
        + a2xx * (b2xx + 2.0 * b3xx * y)
        - 2.0 * a3xy * (b1xy - b3xy * x**2 + b2xy * y + b3xy * y**2)
        + 2.0
        * (
            x * (-(a3yy * b2yy) + a2xy * b3xy - 8.0 * at * (-1.0 + x) ** 2 * x)
            - b3yy * (a2yy + 2.0 * a3yy * x) * y
            + 8.0 * at * (-1.0 + 6.0 * (x - 2.0 * x**3 + x**4)) * y**2
            + 8.0 * at * (1.0 + 6.0 * (-1.0 + x) * x) * y**4
            + a3xx * x * (b2xx + 2.0 * b3xx * y)
        )
    )

    return (
        betann * exp_at * betann_term + Rey * (convective_term + exp_at * exp_term)
    ) / (exp_2at * Rey)


def evaluate_case_02_tftxx(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxx`` source term for MMS Case 02."""
    c = _case02_common(grid, t, parameters)
    p: Case02SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    Axx: FloatArray = c["Axx"]  # type: ignore[assignment]
    Bxx: FloatArray = c["Bxx"]  # type: ignore[assignment]
    Axy: FloatArray = c["Axy"]  # type: ignore[assignment]
    Bxy: FloatArray = c["Bxy"]  # type: ignore[assignment]
    tensor_sum: FloatArray = c["tensor_sum"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    epsylon = p.epsilon
    xi = p.xi
    alphaG = p.alpha_g
    at = p.at

    velocity_term = (
        64.0 * (-1.0 + betann) * exp_at * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y * (-1.0 + 2.0 * y**2)
        / Rey
    )

    giesekus_term = alphaG * (1.0 - betann) * Rey * Wi * (Axx**2 * Bxx**2 + Axy**2 * Bxy**2)

    wi_inner = (
        16.0 * (-1.0 + betann) * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Axx * y**2 * (p.b2xx + 2.0 * p.b3xx * y)
        * (-1.0 + y**2)
        + at * (-1.0 + betann) * exp_at * Axx * Bxx
        - 16.0 * (-1.0 + betann) * (-1.0 + x) ** 2 * x**2 * (p.a2xx + 2.0 * p.a3xx * x) * y * (-1.0 + 2.0 * y**2)
        * Bxx
        + 64.0 * (-1.0 + betann) * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Axx * y * (-1.0 + 2.0 * y**2) * Bxx
        + 32.0 * (-1.0 + betann) * (-1.0 + x) ** 2 * x**2 * Axy * (-1.0 + 6.0 * y**2) * Bxy
        + 16.0
        * (1.0 - betann)
        * xi
        * (
            4.0 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Axx * y * (-1.0 + 2.0 * y**2) * Bxx
            - (1.0 + 6.0 * (-1.0 + x) * x) * Axy * y**2 * (-1.0 + y**2) * Bxy
            + (-1.0 + x) ** 2 * x**2 * Axy * (-1.0 + 6.0 * y**2) * Bxy
        )
    )

    relaxation_term = (1.0 - betann) * Axx * Bxx * (exp_at + epsylon * Rey * Wi * tensor_sum)

    return (velocity_term + giesekus_term + Wi * wi_inner + relaxation_term) / exp_2at


def evaluate_case_02_tftxy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTxy`` source term for MMS Case 02."""
    c = _case02_common(grid, t, parameters)
    p: Case02SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    Axx: FloatArray = c["Axx"]  # type: ignore[assignment]
    Bxx: FloatArray = c["Bxx"]  # type: ignore[assignment]
    Axy: FloatArray = c["Axy"]  # type: ignore[assignment]
    Bxy: FloatArray = c["Bxy"]  # type: ignore[assignment]
    Ayy: FloatArray = c["Ayy"]  # type: ignore[assignment]
    Byy: FloatArray = c["Byy"]  # type: ignore[assignment]
    tensor_sum: FloatArray = c["tensor_sum"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    xi = p.xi
    alphaG = p.alpha_g
    epsylon = p.epsilon
    at = p.at

    velocity_term = (
        -16.0
        * (-1.0 + betann)
        * exp_at
        * (
            (-1.0 + x) ** 2 * x**2
            + (-1.0 - 6.0 * (-1.0 + x) * x * (1.0 + (-1.0 + x) * x)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * y**4
        )
        / Rey
    )

    giesekus_term = alphaG * (1.0 - betann) * Rey * Wi * Axy * Bxy * (Axx * Bxx + Ayy * Byy)

    relaxation_term = -(-1.0 + betann) * Axy * Bxy * (exp_at + epsylon * Rey * Wi * tensor_sum)

    wi_inner = (
        16.0 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Axy * y**2 * (p.b2xy + 2.0 * p.b3xy * y) * (-1.0 + y**2)
        - 16.0 * (1.0 + 6.0 * (-1.0 + x) * x) * Axx * y**2 * (-1.0 + y**2) * Bxx
        + at * exp_at * Axy * Bxy
        - 16.0 * (-1.0 + x) ** 2 * x**2 * (p.a2xy + 2.0 * p.a3xy * x) * y * (-1.0 + 2.0 * y**2) * Bxy
        + 16.0 * (-1.0 + x) ** 2 * x**2 * Ayy * (-1.0 + 6.0 * y**2) * Byy
        + 8.0
        * xi
        * (
            (-1.0 + x) ** 2 * x**2
            + (-1.0 - 6.0 * (-1.0 + x) * x * (1.0 + (-1.0 + x) * x)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * y**4
        )
        * tensor_sum
    )

    return (
        velocity_term + giesekus_term + relaxation_term + (-1.0 + betann) * Wi * wi_inner
    ) / exp_2at


def evaluate_case_02_tftyy(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfTyy`` source term for MMS Case 02."""
    c = _case02_common(grid, t, parameters)
    p: Case02SourceTermParameters = c["p"]  # type: ignore[assignment]
    x: FloatArray = c["x"]  # type: ignore[assignment]
    y: FloatArray = c["y"]  # type: ignore[assignment]
    exp_at: FloatArray = c["exp_at"]  # type: ignore[assignment]
    exp_2at: FloatArray = c["exp_2at"]  # type: ignore[assignment]
    Axy: FloatArray = c["Axy"]  # type: ignore[assignment]
    Bxy: FloatArray = c["Bxy"]  # type: ignore[assignment]
    Ayy: FloatArray = c["Ayy"]  # type: ignore[assignment]
    Byy: FloatArray = c["Byy"]  # type: ignore[assignment]
    tensor_sum: FloatArray = c["tensor_sum"]  # type: ignore[assignment]

    Rey = p.reynolds
    Wi = p.wi
    betann = p.beta
    xi = p.xi
    alphaG = p.alpha_g
    epsylon = p.epsilon
    at = p.at

    velocity_term = (
        -64.0 * (-1.0 + betann) * exp_at * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y * (-1.0 + 2.0 * y**2) / Rey
    )

    giesekus_term = alphaG * (1.0 - betann) * Rey * Wi * (Axy**2 * Bxy**2 + Ayy**2 * Byy**2)

    relaxation_term = (1.0 - betann) * Ayy * Byy * (exp_at + epsylon * Rey * Wi * tensor_sum)

    wi_inner = (
        16.0 * (-1.0 + betann) * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Ayy * y**2 * (p.b2yy + 2.0 * p.b3yy * y)
        * (-1.0 + y**2)
        - 32.0 * (-1.0 + betann) * (1.0 + 6.0 * (-1.0 + x) * x) * Axy * y**2 * (-1.0 + y**2) * Bxy
        + at * (-1.0 + betann) * exp_at * Ayy * Byy
        - 16.0 * (-1.0 + betann) * (-1.0 + x) ** 2 * x**2 * (p.a2yy + 2.0 * p.a3yy * x) * y * (-1.0 + 2.0 * y**2)
        * Byy
        - 64.0 * (-1.0 + betann) * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Ayy * y * (-1.0 + 2.0 * y**2) * Byy
        + 16.0
        * (1.0 - betann)
        * xi
        * (
            -((1.0 + 6.0 * (-1.0 + x) * x) * Axy * y**2 * (-1.0 + y**2) * Bxy)
            + (-1.0 + x) ** 2 * x**2 * Axy * (-1.0 + 6.0 * y**2) * Bxy
            - 4.0 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * Ayy * y * (-1.0 + 2.0 * y**2) * Byy
        )
    )

    return (velocity_term + giesekus_term + Wi * wi_inner + relaxation_term) / exp_2at


def evaluate_case_02_source_terms(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02SourceTermParameters | None = None,
) -> Case02SourceTerms:
    """Evaluate all Case 02 source terms."""
    return Case02SourceTerms(
        tfwz=evaluate_case_02_tfwz(grid, t=t, parameters=parameters),
        tfTxx=evaluate_case_02_tftxx(grid, t=t, parameters=parameters),
        tfTxy=evaluate_case_02_tftxy(grid, t=t, parameters=parameters),
        tfTyy=evaluate_case_02_tftyy(grid, t=t, parameters=parameters),
    )


__all__ = [
    "CASE02_FORTRAN_ARRAYS",
    "Case02SourceTermParameters",
    "Case02SourceTerms",
    "evaluate_case_02_source_terms",
    "evaluate_case_02_tftxx",
    "evaluate_case_02_tftxy",
    "evaluate_case_02_tftyy",
    "evaluate_case_02_tfwz",
]
