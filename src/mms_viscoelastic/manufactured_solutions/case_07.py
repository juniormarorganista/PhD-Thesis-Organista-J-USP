"""Manufactured solution: Fortran ``case_07.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_07(t)``. Unlike every other ported case,
the velocity/vorticity/stream-function fields here decay in time through a
hyperbolic tangent, ``-1 + tanh(4 - 8 * at * t)``, instead of the
``exp(-at*t)`` decay used elsewhere. The polymeric stress components are
spatially uniform (constant in ``x``/``y``, only decaying via ``exp(-at*t)``),
the same ansatz as case 01 (fixed targets ``c1``, ``c2``, ``c3`` from
``par.nn``, scaled by ``(1 - beta)``).

``boundaries_condictions_case_07`` re-evaluates this same exact solution
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
from mms_viscoelastic.parameters import FlowParameters

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case07Parameters:
    """Parameters specific to the exact MMS case 07 fields.

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
    ) -> "Case07Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, c1=c1, c2=c2, c3=c3)

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def evaluate_case_07_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case07Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_07.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-07 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case07Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    tanh_term = -1.0 + np.tanh(4.0 - 8.0 * params.at * t)

    x_poly = (-1.0 + x) ** 2 * x**2

    u = -8.0 * x_poly * y**2 * tanh_term

    v = 16.0 * (-1.0 + x) * x * (-1.0 + 2.0 * x) * y**3 * tanh_term / 3.0

    wz = -16.0 * y * (y**2 + 3.0 * (-1.0 + x) * x * ((-1.0 + x) * x + 2.0 * y**2)) * tanh_term / 3.0

    psi = -8.0 * x_poly * y**3 * tanh_term / 3.0

    # txx/txy/tyy are spatially uniform in case 07 (constant targets c1/c2/c3,
    # only decaying in time via exp(-at*t), same ansatz as case 01); broadcast
    # to the grid shape so all Field2D components share the same shape, as
    # required by Field2D.as_feature_matrix.
    txx = np.full_like(x, -((-1.0 + params.beta) * params.c1) / decay)
    txy = np.full_like(x, -((-1.0 + params.beta) * params.c2) / decay)
    tyy = np.full_like(x, -((-1.0 + params.beta) * params.c3) / decay)

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_07_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 07."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
