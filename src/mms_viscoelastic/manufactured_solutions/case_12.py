"""Manufactured solution: Fortran ``case_12.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_12(t)``.  It is intentionally focused on the
exact solution first; the very long source terms from ``source_term_case_12(t)``
should be ported in a later phase after this field-level translation is tested.

Array convention
----------------
Fortran stores fields as ``field(i, j)``.  NumPy stores the same physical field
as arrays with shape ``(ny, nx)``, where axis 0 is ``y`` and axis 1 is ``x``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.parameters import FlowParameters

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case12Coefficients:
    """Trigonometric stress coefficients used by MMS case 12.

    The defaults reproduce the values defined in the original ``par.nn`` file:
    all ``a`` and ``b`` coefficients are equal to one.
    """

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
class Case12Parameters:
    """Parameters specific to the exact MMS case 12 fields."""

    beta: float = 0.5
    at: float = 0.0
    coefficients: Case12Coefficients = Case12Coefficients()

    @classmethod
    def from_flow_parameters(
        cls,
        flow: FlowParameters,
        *,
        at: float = 0.0,
        coefficients: Case12Coefficients | None = None,
    ) -> "Case12Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, coefficients=coefficients or Case12Coefficients())

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def _stress_factor_x(x: FloatArray, a1: float, a2: float, a3: float) -> FloatArray:
    return a1 + a3 * np.cos(np.pi * x) + a2 * np.sin(np.pi * x)


def _stress_factor_y(y: FloatArray, b1: float, b2: float, b3: float) -> FloatArray:
    return b1 + b3 * np.cos(np.pi * y) + b2 * np.sin(np.pi * y)


def evaluate_case_12_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_12.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid.  The returned arrays have shape ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-12 parameters.  If omitted, defaults mirror the thesis code for
        ``beta`` and trigonometric stress coefficients.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and ``tyy``.
    """
    params = parameters or Case12Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    pi = np.pi
    decay = np.exp(params.at * t)
    c = params.coefficients

    sin_px = np.sin(pi * x)
    cos_px = np.cos(pi * x)
    sin_py = np.sin(pi * y)
    cos_py = np.cos(pi * y)

    x_poly = (-1.0 + x) * x

    u = x_poly * sin_px * (pi * y * cos_py + sin_py) / decay

    v = -(
        y
        * (pi * x_poly * cos_px + (-1.0 + 2.0 * x) * sin_px)
        * sin_py
    ) / decay

    wz = (
        2.0
        * (
            pi * x_poly * cos_py * sin_px
            + y
            * (
                pi * (-1.0 + 2.0 * x) * cos_px
                + (1.0 - pi**2 * x_poly) * sin_px
            )
            * sin_py
        )
    ) / decay

    psi = x_poly * y * sin_px * sin_py / decay

    txx = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1xx, c.a2xx, c.a3xx)
        * _stress_factor_y(y, c.b1xx, c.b2xx, c.b3xx)
        / decay
    )
    txy = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1xy, c.a2xy, c.a3xy)
        * _stress_factor_y(y, c.b1xy, c.b2xy, c.b3xy)
        / decay
    )
    tyy = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1yy, c.a2yy, c.a3yy)
        * _stress_factor_y(y, c.b1yy, c.b2yy, c.b3yy)
        / decay
    )

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_12_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 12."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
