"""Manufactured solution: Fortran ``case_08.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_08(t)``. Like case 07, the velocity,
vorticity and stream-function fields decay in time through a hyperbolic
tangent, ``1 - tanh(4 - 8 * at * t)``, instead of ``exp(-at*t)``, but with
a different (purely polynomial) spatial structure. The polymeric stress
components are quadratic polynomials in ``x`` and ``y``, the same ansatz
as case 02, reusing
:class:`~mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients`.

``boundaries_condictions_case_08`` re-evaluates this same exact solution
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
class Case08Parameters:
    """Parameters specific to the exact MMS case 08 fields."""

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
    ) -> "Case08Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, coefficients=coefficients or Case12Coefficients())

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def _quadratic_x(x: FloatArray, a1: float, a2: float, a3: float) -> FloatArray:
    return a1 + a2 * x + a3 * x**2


def _quadratic_y(y: FloatArray, b1: float, b2: float, b3: float) -> FloatArray:
    return b1 + b2 * y + b3 * y**2


def evaluate_case_08_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case08Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_08.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-08 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case08Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = 1.0 - np.tanh(4.0 - 8.0 * params.at * t)
    x_poly = (1.0 - x) ** 2 * x**2

    u = 8.0 * x_poly * y**2 * decay

    v = (
        -16.0 * (1.0 - x) ** 2 * x * y**3 * decay / 3.0
        + 16.0 * (1.0 - x) * x**2 * y**3 * decay / 3.0
    )

    wz = (
        16.0 * x_poly * y * decay
        + 16.0 * (1.0 - x) ** 2 * y**3 * decay / 3.0
        - 64.0 * (1.0 - x) * x * y**3 * decay / 3.0
        + 16.0 * x**2 * y**3 * decay / 3.0
    )

    psi = -16.0 * (0.5 * x**2 - x**3 + 0.5 * x**4) * y**3 * (-1.0 + np.tanh(4.0 - 8.0 * params.at * t)) / 3.0

    c = params.coefficients
    decay_exp = np.exp(params.at * t)
    txx = (1.0 - params.beta) * _quadratic_x(x, c.a1xx, c.a2xx, c.a3xx) * _quadratic_y(y, c.b1xx, c.b2xx, c.b3xx) / decay_exp
    txy = (1.0 - params.beta) * _quadratic_x(x, c.a1xy, c.a2xy, c.a3xy) * _quadratic_y(y, c.b1xy, c.b2xy, c.b3xy) / decay_exp
    tyy = (1.0 - params.beta) * _quadratic_x(x, c.a1yy, c.a2yy, c.a3yy) * _quadratic_y(y, c.b1yy, c.b2yy, c.b3yy) / decay_exp

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_08_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 08."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
