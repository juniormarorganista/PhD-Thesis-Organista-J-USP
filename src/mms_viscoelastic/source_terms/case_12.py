"""Source terms for MMS Case 12.

This module starts the algebraic Python port of the Fortran routine:

``src/fortran_original/4Order_LPOG_Euler/case_12.f::source_term_case_12(t)``

Current implementation status:

- ``tfwz``: ported to NumPy;
- ``tfTxx``: pending;
- ``tfTxy``: pending;
- ``tfTyy``: pending.

The full source-term evaluator remains protected by default. It raises
``NotImplementedError`` unless explicitly called with ``allow_partial=True``.
This prevents incomplete source terms from being used silently in verification
or machine-learning workflows.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]

CASE12_FORTRAN_ARRAYS: tuple[str, ...] = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
CASE12_PORTED_ARRAYS: tuple[str, ...] = ("tfwz",)
CASE12_PENDING_ARRAYS: tuple[str, ...] = ("tfTxx", "tfTxy", "tfTyy")


@dataclass(frozen=True)
class Case12SourceTermParameters:
    """Physical and constitutive parameters required by Case 12 source terms.

    Defaults are intentionally simple and mirror the early Python case-12
    defaults: unit trigonometric coefficients, ``Re = 1`` and ``beta = 0.5``.
    """

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
class Case12SourceTerms:
    """Container for MMS Case 12 source-term arrays.

    Pending components are represented by ``None`` until their full algebraic
    Fortran expressions are ported and validated.
    """

    tfwz: Optional[FloatArray] = None
    tfTxx: Optional[FloatArray] = None
    tfTxy: Optional[FloatArray] = None
    tfTyy: Optional[FloatArray] = None

    def available_components(self) -> tuple[str, ...]:
        """Return source-term components currently stored as arrays."""

        names = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
        return tuple(name for name in names if getattr(self, name) is not None)

    def pending_components(self) -> tuple[str, ...]:
        """Return source-term components not yet stored as arrays."""

        names = ("tfwz", "tfTxx", "tfTxy", "tfTyy")
        return tuple(name for name in names if getattr(self, name) is None)


def available_case12_source_terms() -> tuple[str, ...]:
    """Return the Fortran source-term arrays targeted by the Python port."""

    return CASE12_FORTRAN_ARRAYS


def ported_case12_source_terms() -> tuple[str, ...]:
    """Return source-term arrays already ported to Python."""

    return CASE12_PORTED_ARRAYS


def pending_case12_source_terms() -> tuple[str, ...]:
    """Return source-term arrays still pending."""

    return CASE12_PENDING_ARRAYS


def case12_source_term_status() -> dict[str, str]:
    """Return current implementation status for each source-term component."""

    return {
        "tfwz": "ported to NumPy; regression against Fortran reference pending",
        "tfTxx": "mapped from Fortran, algebraic Python port pending",
        "tfTxy": "mapped from Fortran, algebraic Python port pending",
        "tfTyy": "mapped from Fortran, algebraic Python port pending",
    }


def evaluate_case_12_tfwz(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
) -> FloatArray:
    """Evaluate the ported ``tfwz`` source term for MMS Case 12.

    The expression is a direct NumPy translation of the Fortran assignment:

    ``tfwz(i,j) = ...``

    in ``source_term_case_12(t)``.

    Notes
    -----
    Fortran stores ``tfwz(i,j)``. Python returns an array with shape
    ``(grid.ny, grid.nx)``, so the same physical value is accessed as
    ``tfwz[j, i]`` using zero-based Python indices.
    """

    params = parameters or Case12SourceTermParameters()
    params.validate()

    x = grid.X
    y = grid.Y
    pi = np.pi

    Rey = params.reynolds
    betann = params.beta
    at = params.at

    a1xx = params.a1xx
    a2xx = params.a2xx
    a3xx = params.a3xx
    b1xx = params.b1xx
    b2xx = params.b2xx
    b3xx = params.b3xx

    a1xy = params.a1xy
    a2xy = params.a2xy
    a3xy = params.a3xy
    b1xy = params.b1xy
    b2xy = params.b2xy
    b3xy = params.b3xy

    a1yy = params.a1yy
    a2yy = params.a2yy
    a3yy = params.a3yy
    b1yy = params.b1yy
    b2yy = params.b2yy
    b3yy = params.b3yy

    exp_at = np.exp(at * t)
    exp_2at = np.exp(2.0 * at * t)

    cos_px = np.cos(pi * x)
    sin_px = np.sin(pi * x)
    cos_py = np.cos(pi * y)
    sin_py = np.sin(pi * y)

    cos_2px = np.cos(2.0 * pi * x)
    sin_2px = np.sin(2.0 * pi * x)
    cos_2py = np.cos(2.0 * pi * y)
    sin_2py = np.sin(2.0 * pi * y)

    x_poly = (-1.0 + x) * x

    term_base = (
        -(Rey * y)
        + 2.0 * Rey * x * y
        + Rey * y * cos_2px
        - 2.0 * Rey * x * y * cos_2px
        + 2.0 * pi**2 * Rey * x * y * cos_2px
        - 6.0 * pi**2 * Rey * x**2 * y * cos_2px
        + 4.0 * pi**2 * Rey * x**3 * y * cos_2px
        + pi * Rey * y * sin_2px
        - 2.0 * pi * Rey * x * y * sin_2px
        + 2.0 * pi * Rey * x**2 * y * sin_2px
        - 2.0 * pi**3 * Rey * x**2 * y * sin_2px
        + 4.0 * pi**3 * Rey * x**3 * y * sin_2px
        - 2.0 * pi**3 * Rey * x**4 * y * sin_2px
    )

    term_cos2y = (
        Rey
        * y
        * cos_2py
        * (
            -(
                (-1.0 + 2.0 * x)
                * (1.0 + 2.0 * pi**2 * x_poly - cos_2px)
            )
            + pi
            * (-1.0 - 2.0 * x_poly)
            * sin_2px
        )
    )

    term_cosx = (
        2.0
        * exp_at
        * pi
        * cos_px
        * (
            -(a3xy * b1xy * (-1.0 + betann) * pi * Rey)
            + pi
            * (
                (a2xx * b2xx - a2yy * b2yy) * Rey
                + betann
                * (
                    -8.0
                    - a2xx * b2xx * Rey
                    + a2yy * b2yy * Rey
                    + 16.0 * x
                )
            )
            * cos_py
            + (
                (a2xx * b3xx - a2yy * b3yy)
                * (-1.0 + betann)
                * pi
                * Rey
                - 2.0
                * (4.0 * betann * pi**2 - at * Rey)
                * (-1.0 + 2.0 * x)
                * y
            )
            * sin_py
        )
    )

    term_exp = exp_at * (
        -2.0
        * a2xy
        * b1xy
        * (-1.0 + betann)
        * pi**2
        * Rey
        * sin_px
        + 2.0
        * pi
        * cos_py
        * (
            a1xy * b3xy * (-1.0 + betann) * pi * Rey
            + (
                Rey
                * (
                    -(a3xx * b2xx * pi)
                    + a3yy * b2yy * pi
                    + 2.0 * at * x_poly
                )
                + betann
                * (
                    8.0
                    + a3xx * b2xx * pi * Rey
                    - pi
                    * (
                        a3yy * b2yy * Rey
                        + 8.0 * pi * x_poly
                    )
                )
            )
            * sin_px
        )
        + 2.0
        * (
            a1xy * b2xy * (-1.0 + betann) * pi**2 * Rey
            + (
                -(
                    (a3xx * b3xx - a3yy * b3yy)
                    * (-1.0 + betann)
                    * pi**2
                    * Rey
                )
                + 2.0
                * (
                    at * Rey
                    + pi**2
                    * (
                        -(at * Rey * x_poly)
                        + 2.0 * betann * (-4.0 + pi**2 * x_poly)
                    )
                )
                * y
            )
            * sin_px
        )
        * sin_py
    )

    term_sin2y = (
        pi
        * Rey
        * (
            (-1.0 + 2.0 * x)
            * (y**2 + x_poly * (-1.0 + 2.0 * pi**2 * y**2))
            + (-1.0 + 2.0 * x)
            * (x_poly - y**2)
            * cos_2px
            - pi
            * (
                ((-1.0 + x) ** 2) * x**2
                + (-1.0 - 2.0 * x_poly) * y**2
            )
            * sin_2px
        )
        * sin_2py
    )

    numerator = term_base + term_cos2y + term_cosx + term_exp + term_sin2y

    return numerator / (2.0 * exp_2at * Rey)


def evaluate_case_12_source_terms(
    grid: Grid2D | None = None,
    t: float = 0.0,
    parameters: Case12SourceTermParameters | None = None,
    *,
    allow_partial: bool = False,
) -> Case12SourceTerms:
    """Evaluate MMS Case 12 source terms.

    By default this function still raises ``NotImplementedError`` because only
    ``tfwz`` has been ported. Use ``allow_partial=True`` when the caller
    explicitly accepts a partial result with pending tensor terms set to
    ``None``.
    """

    if not allow_partial:
        raise NotImplementedError(
            "MMS Case 12 full source terms are not fully ported yet. "
            "Only tfwz is currently available through evaluate_case_12_tfwz(...). "
            "Use evaluate_case_12_source_terms(..., allow_partial=True) only "
            "when a partial source-term object is intentional."
        )

    if grid is None:
        raise ValueError("grid is required when allow_partial=True.")

    tfwz = evaluate_case_12_tfwz(grid, t=t, parameters=parameters)
    return Case12SourceTerms(tfwz=tfwz, tfTxx=None, tfTxy=None, tfTyy=None)


__all__ = [
    "CASE12_FORTRAN_ARRAYS",
    "CASE12_PENDING_ARRAYS",
    "CASE12_PORTED_ARRAYS",
    "Case12SourceTermParameters",
    "Case12SourceTerms",
    "available_case12_source_terms",
    "case12_source_term_status",
    "evaluate_case_12_source_terms",
    "evaluate_case_12_tfwz",
    "pending_case12_source_terms",
    "ported_case12_source_terms",
]
