"""Manufactured solution: Fortran ``case_01.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_01(t)``. Unlike case 12, the polymeric
tensor components of case 01 are spatially uniform (constant in ``x``/``y``,
only decaying in time): they reproduce the fixed targets ``c1=1``, ``c2=2``,
``c3=3`` from ``par.nn``, scaled by ``(1 - beta)``.

``boundaries_condictions_case_01`` in the original Fortran re-evaluates this
same exact solution restricted to the boundary nodes; it has no separate
formula and is therefore not ported as a distinct function, following the
same convention already used for case 12 (see
``mms_viscoelastic.numerics.time_marching``).

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
class Case01Parameters:
    """Parameters specific to the exact MMS case 01 fields.

    ``c1``, ``c2`` and ``c3`` are the fixed target stress values from
    ``par.nn`` (``parameter (c1 = 1.0d0)``, etc.), shared across the MMS
    cases that use a spatially uniform stress field.
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
    ) -> "Case01Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, c1=c1, c2=c2, c3=c3)

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def evaluate_case_01_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case01Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_01.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-01 parameters. If omitted, defaults mirror ``par.nn``.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and
        ``tyy``.
    """
    params = parameters or Case01Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    decay = np.exp(params.at * t)

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

    # txx/txy/tyy are spatially uniform in case 01 (constant targets c1/c2/c3,
    # only decaying in time); broadcast to the grid shape so all Field2D
    # components share the same shape, as required by Field2D.as_feature_matrix.
    txx = np.full_like(x, -((-1.0 + params.beta) * params.c1) / decay)
    txy = np.full_like(x, -((-1.0 + params.beta) * params.c2) / decay)
    tyy = np.full_like(x, -((-1.0 + params.beta) * params.c3) / decay)

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_01_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 01."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
