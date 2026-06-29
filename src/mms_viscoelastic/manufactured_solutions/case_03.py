"""Manufactured solution: Fortran ``case_03.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_03(t)``. The velocity/vorticity/stream
function fields (``u``, ``v``, ``wz``, ``psi``) are algebraically identical
to cases 01 and 02. The polymeric stress components use the same
trigonometric ansatz as case 12
(``(1 - beta) * (a1 + a3 cos(pi x) + a2 sin(pi x)) * (b1 + b3 cos(pi y) +
b2 sin(pi y))``), reusing
:class:`~mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients`
since the field names and defaults (all coefficients equal to one) match
``par.nn`` exactly.

``boundaries_condictions_case_03`` re-evaluates this same exact solution
restricted to the boundary nodes; it has no separate formula and is
therefore not ported as a distinct function (same convention as cases 01
and 02).

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
class Case03Parameters:
    """Parameters specific to the exact MMS case 03 fields."""

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
    ) -> "Case03Parameters":
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


def evaluate_case_03_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case03Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_03.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-03 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case03Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    c = params.coefficients

    x_quartic = x**2 - 2.0 * x**3 + x**4
    x_quad = 1.0 - 3.0 * x + 2.0 * x**2
    y_quartic = -0.5 * y**2 + 0.5 * y**4

    u = 8.0 * x_quartic * (-2.0 * y + 4.0 * y**3) / decay

    v = -32.0 * x * x_quad * y_quartic / decay

    wz = (
        8.0 * x_quartic * (-2.0 + 12.0 * y**2) / decay
        + 32.0 * x * (-3.0 + 4.0 * x) * y_quartic / decay
        + 32.0 * x_quad * y_quartic / decay
    )

    psi = 16.0 * (0.5 * x**2 - x**3 + 0.5 * x**4) * y**2 * (-1.0 + y**2) / decay

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


def case_03_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 03."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
