"""Manufactured solution: Fortran ``case_10.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_10(t)``. The velocity, vorticity and
stream-function fields (``u``, ``v``, ``wz``, ``psi``) follow the same
trigonometric-times-polynomial ansatz as case 12
(``u = x_poly * sin(pi x) * (pi y cos(pi y) + sin(pi y)) / decay``, with
``x_poly = (-1 + x) * x``). The polymeric stress components are
spatially uniform, the same ansatz as case 01 (fixed targets ``c1``,
``c2``, ``c3`` from ``par.nn``, scaled by ``(1 - beta)``).

``boundaries_condictions_case_10`` additionally calls an external
``boundarie_vorticity_cavity_flow`` subroutine for a one-sided
finite-difference wall-vorticity boundary condition (lines that are
commented out in the original source show the explicit formula); since
the interior exact fields/source terms do not depend on this boundary
treatment, and this module only ports the exact-solution formulas (not
boundary conditions, following the same convention as the other ported
cases), this detail does not affect the port.

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
from mms_viscoelastic.parameters import FlowParameters

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case10Parameters:
    """Parameters specific to the exact MMS case 10 fields.

    ``c1``, ``c2`` and ``c3`` are the fixed target stress values from
    ``par.nn`` (``parameter (c1 = 1.0d0)``, etc.), the same convention
    used by case 01.
    """

    beta: float = 0.5
    at: float = 0.0
    c1: float = 1.0
    c2: float = 2.0
    c3: float = 3.0

    @classmethod
    def from_flow_parameters(
        cls,
        flow: FlowParameters,
        *,
        at: float = 0.0,
        c1: float = 1.0,
        c2: float = 2.0,
        c3: float = 3.0,
    ) -> "Case10Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, c1=c1, c2=c2, c3=c3)

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def evaluate_case_10_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case10Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_10.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-10 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case10Parameters()
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

    # txx/txy/tyy are spatially uniform in case 10 (constant targets
    # c1/c2/c3, only decaying in time), same ansatz as case 01; broadcast
    # to the grid shape so all Field2D components share the same shape, as
    # required by Field2D.as_feature_matrix.
    txx = np.full_like(x, -((-1.0 + params.beta) * params.c1) / decay)
    txy = np.full_like(x, -((-1.0 + params.beta) * params.c2) / decay)
    tyy = np.full_like(x, -((-1.0 + params.beta) * params.c3) / decay)

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_10_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 10."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
