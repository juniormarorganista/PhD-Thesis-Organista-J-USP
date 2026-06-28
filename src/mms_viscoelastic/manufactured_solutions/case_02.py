"""Manufactured solution: Fortran ``case_02.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_02(t)``. The velocity/vorticity/stream
function fields (``u``, ``v``, ``wz``, ``psi``) are algebraically identical
to case 01. The polymeric stress components are quadratic polynomials in
``x`` and ``y`` (``a1 + a2 x + a3 x^2`` times ``b1 + b2 y + b3 y^2``),
reusing :class:`~mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients`
since the field names and defaults (all coefficients equal to one) match
``par.nn`` exactly.

``boundaries_condictions_case_02`` re-evaluates this same exact solution
restricted to the boundary nodes; it has no separate formula and is
therefore not ported as a distinct function (same convention as case 01).

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
class Case02Parameters:
    """Parameters specific to the exact MMS case 02 fields."""

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
    ) -> "Case02Parameters":
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


def evaluate_case_02_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case02Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_02.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-02 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case02Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    c = params.coefficients

    x_poly = (-1.0 + x) ** 2 * x**2

    u = 16.0 * x_poly * y * (-1.0 + 2.0 * y**2) / decay

    v = -16.0 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y**2 * (-1.0 + y**2) / decay

    wz = (
        16.0
        * (
            -x_poly
            + (-1.0 + 6.0 * (x - 2.0 * x**3 + x**4)) * y**2
            + (1.0 + 6.0 * (-1.0 + x) * x) * y**4
        )
        / decay
    )

    psi = 8.0 * x_poly * y**2 * (-1.0 + y**2) / decay

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


def case_02_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 02."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
