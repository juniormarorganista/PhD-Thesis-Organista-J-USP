"""Manufactured solution: Fortran ``case_04.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_04(t)``. The velocity/vorticity/stream
function fields (``u``, ``v``, ``wz``, ``psi``) follow a Taylor-Green-vortex
trigonometric ansatz, distinct from cases 01/02/03 and case 12. The
polymeric stress components are spatially uniform (constant in ``x``/``y``,
only decaying in time), algebraically identical to case 01's ansatz
(fixed targets ``c1``, ``c2``, ``c3`` from ``par.nn``, scaled by
``(1 - beta)``).

``boundaries_condictions_case_04`` re-evaluates this same exact solution
restricted to the boundary nodes; it has no separate formula and is
therefore not ported as a distinct function (same convention as cases 01,
02 and 03).

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
class Case04Parameters:
    """Parameters specific to the exact MMS case 04 fields.

    ``c1``, ``c2`` and ``c3`` are the fixed target stress values from
    ``par.nn`` (``parameter (c1 = 1.0d0)``, etc.), shared across the MMS
    cases that use a spatially uniform stress field (cases 01 and 04).
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
    ) -> "Case04Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, c1=c1, c2=c2, c3=c3)

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def evaluate_case_04_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case04Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_04.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-04 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case04Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)
    pi = np.pi

    cos_px = np.cos(pi * x)
    sin_px = np.sin(pi * x)
    cos_py = np.cos(pi * y)
    sin_py = np.sin(pi * y)

    u = 8.0 * cos_px * sin_py / decay
    v = -8.0 * cos_py * sin_px / decay
    wz = 16.0 * pi * cos_px * cos_py / decay
    psi = -8.0 * cos_px * cos_py / (decay * pi)

    # txx/txy/tyy are spatially uniform in case 04 (constant targets
    # c1/c2/c3, only decaying in time, same ansatz as case 01); broadcast
    # to the grid shape so all Field2D components share the same shape, as
    # required by Field2D.as_feature_matrix.
    txx = np.full_like(x, -((-1.0 + params.beta) * params.c1) / decay)
    txy = np.full_like(x, -((-1.0 + params.beta) * params.c2) / decay)
    tyy = np.full_like(x, -((-1.0 + params.beta) * params.c3) / decay)

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_04_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 04."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
