"""Manufactured solution: Fortran ``case_11.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_11(t)``. The velocity, vorticity and
stream-function fields (``u``, ``v``, ``wz``, ``psi``) are algebraically
identical to case 10 (and case 12). The polymeric stress components are
quadratic polynomials in ``x`` and ``y``, the same ansatz as case 02,
reusing
:class:`~mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients`.

``boundaries_condictions_case_11`` re-evaluates this same exact solution
restricted to the boundary nodes; it has no separate formula and is
therefore not ported as a distinct function (same convention as the
other ported cases).

Array convention
----------------
Fortran stores fields as ``field(i, j)``. NumPy stores the same physical
field as arrays with shape ``(ny, nx)``, where axis 0 is ``y`` and axis 1 is
``x``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions.case_12 import Case12Coefficients
from mms_viscoelastic.parameters import FlowParameters

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case11Parameters:
    """Parameters specific to the exact MMS case 11 fields."""

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
    ) -> "Case11Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, coefficients=coefficients or Case12Coefficients())

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def _quadratic_x(x: FloatArray, a1: float, a2: float, a3: float) -> FloatArray:
    return a1 + x * (a2 + a3 * x)


def _quadratic_y(y: FloatArray, b1: float, b2: float, b3: float) -> FloatArray:
    return b1 + y * (b2 + b3 * y)


def evaluate_case_11_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case11Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_11.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-11 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case11Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    pi = np.pi

    cos_px = np.cos(pi * x)
    sin_px = np.sin(pi * x)
    cos_py = np.cos(pi * y)
    sin_py = np.sin(pi * y)
    x_poly = (-1.0 + x) * x

    u = x_poly * sin_px * (pi * y * cos_py + sin_py) / decay

    v = -(y * (pi * x_poly * cos_px + (-1.0 + 2.0 * x) * sin_px) * sin_py) / decay

    wz = (
        2.0
        * (
            pi * x_poly * cos_py * sin_px
            + y * (pi * (-1.0 + 2.0 * x) * cos_px + (1.0 - pi**2 * x_poly) * sin_px) * sin_py
        )
        / decay
    )

    psi = x_poly * y * sin_px * sin_py / decay

    c = params.coefficients
    txx = (
        -(-1.0 + params.beta)
        * _quadratic_x(x, c.a1xx, c.a2xx, c.a3xx)
        * _quadratic_y(y, c.b1xx, c.b2xx, c.b3xx)
        / decay
    )
    txy = (
        -(-1.0 + params.beta)
        * _quadratic_x(x, c.a1xy, c.a2xy, c.a3xy)
        * _quadratic_y(y, c.b1xy, c.b2xy, c.b3xy)
        / decay
    )
    tyy = (
        -(-1.0 + params.beta)
        * _quadratic_x(x, c.a1yy, c.a2yy, c.a3yy)
        * _quadratic_y(y, c.b1yy, c.b2yy, c.b3yy)
        / decay
    )

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_11_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 11."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
