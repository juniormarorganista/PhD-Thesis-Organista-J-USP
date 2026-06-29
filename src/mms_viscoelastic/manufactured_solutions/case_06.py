"""Manufactured solution: Fortran ``case_06.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_06(t)``. The velocity/vorticity/stream
function fields (``u``, ``v``, ``wz``, ``psi``) follow the same
Taylor-Green-vortex trigonometric ansatz as cases 04 and 05. The polymeric
stress components use the same trigonometric ansatz as case 03 (and case
12), reusing
:class:`~mms_viscoelastic.manufactured_solutions.case_12.Case12Coefficients`.

``boundaries_condictions_case_06`` re-evaluates this same exact solution
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
class Case06Parameters:
    """Parameters specific to the exact MMS case 06 fields."""

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
    ) -> "Case06Parameters":
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


def evaluate_case_06_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case06Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_06.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-06 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case06Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    pi = np.pi
    c = params.coefficients

    cos_px = np.cos(pi * x)
    sin_px = np.sin(pi * x)
    cos_py = np.cos(pi * y)
    sin_py = np.sin(pi * y)

    u = 8.0 * cos_px * sin_py / decay
    v = -8.0 * cos_py * sin_px / decay
    wz = 16.0 * pi * cos_px * cos_py / decay
    psi = -8.0 * cos_px * cos_py / (decay * pi)

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


def case_06_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 06."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
