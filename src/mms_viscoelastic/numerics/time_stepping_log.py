"""Explicit Euler time stepping for the log-conformation formulation.

Ports the right-hand sides used by ``subroutine drv_Psi`` and
``subroutine nlterms_logPsi`` (``src/fortran_original/4Order_LPOG_Euler/drvs.f``)
and the explicit Euler update used by ``subroutine solver_log_euler``
(``ty_sim == 22`` in ``src/fortran_original/4Order_LPOG_Euler/euler_rk.f``),
for the governing variables ``wz``, ``Psixx``, ``Psixy`` and ``Psiyy``.

As in ``drv_Psi``, the polymeric stress ``(Txx, Txy, Tyy)`` used by the
``wz`` right-hand side is *not* tracked independently here: it is recovered
at every evaluation from the current ``Psi`` state via
:func:`mms_viscoelastic.models.log_conformation.log_conformation_to_stress`
(the Fortran ``call Psi_to_T``).

This path is not discussed in the thesis text and is ported only for
source-code fidelity (Etapa 15;
see ``docs/case-12-log-conformation.md``). As of Etapa 17, the derivatives
use the same 4th-order LPOG compact operators as the non-log formulation
(:mod:`mms_viscoelastic.numerics.time_stepping`), matching the Fortran
``drv_Psi``'s ``derx``/``dery`` calls
(``src/fortran_original/4Order_LPOG_Euler/derivs.f``).

Only interior grid points (``[1:-1, 1:-1]``) are updated; boundary values
are copied unchanged. Re-imposing Dirichlet boundary conditions, the
log-conformation boundary projection (``cc_Psi``) and recomputing
``psi``/``ux``/``uy`` is left to a multi-step time-marching solver, as for
the non-log formulation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.models.log_conformation import (
    kernel_conformation,
    log_conformation_to_stress,
)
from mms_viscoelastic.numerics.derivatives import compact_derivative_x as derivative_x
from mms_viscoelastic.numerics.derivatives import compact_derivative_y as derivative_y
from mms_viscoelastic.numerics.time_stepping import EulerStepParameters, wz_time_derivative

FloatArray = NDArray[np.float64]


def _etap(parameters: EulerStepParameters) -> float:
    """Return ``etap = (1 - beta) / (Rey * Wi)``, matching ``par.nn``."""

    return (1.0 - parameters.beta) / (parameters.reynolds * parameters.wi)


def log_conformation_time_derivatives(
    wz: FloatArray,
    psixx: FloatArray,
    psixy: FloatArray,
    psiyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfwz: FloatArray,
    tf_psixx: FloatArray,
    tf_psixy: FloatArray,
    tf_psiyy: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return ``(dvz, dvPsixx, dvPsixy, dvPsiyy)``.

    Port of ``drv_Psi`` (``drvs.f``): recovers ``(Txx, Txy, Tyy)`` from
    ``Psi`` (``Psi_to_T``), evaluates the rotation rate ``Omegalog`` and the
    tensors ``B``, ``M`` (``kernel_conformation``), and assembles the
    transport right-hand sides for ``wz`` and ``Psi``.
    """

    parameters.validate()
    dx, dy = grid.dx, grid.dy
    etap = _etap(parameters)
    etapinv = 1.0 / etap

    txx, txy, tyy = log_conformation_to_stress(psixx, psixy, psiyy, etap=etap)

    dudx = derivative_x(ux, dx)
    dudy = derivative_y(ux, dy)
    dvdx = derivative_x(uy, dx)
    dvdy = derivative_y(uy, dy)

    omegalog, bxx, bxy, byy, mxx, mxy, myy = kernel_conformation(
        txx, txy, tyy, dudx, dudy, dvdx, dvdy, etapinv=etapinv, alpha_g=parameters.alpha_g
    )

    dvz = wz_time_derivative(wz, txx, txy, tyy, ux, uy, tfwz, grid, parameters)

    u_psixx = ux * psixx
    v_psixx = uy * psixx
    u_psixy = ux * psixy
    v_psixy = uy * psixy
    u_psiyy = ux * psiyy
    v_psiyy = uy * psiyy

    du_psixx_dx = derivative_x(u_psixx, dx)
    dv_psixx_dy = derivative_y(v_psixx, dy)
    du_psixy_dx = derivative_x(u_psixy, dx)
    dv_psixy_dy = derivative_y(v_psixy, dy)
    du_psiyy_dx = derivative_x(u_psiyy, dx)
    dv_psiyy_dy = derivative_y(v_psiyy, dy)

    wi = parameters.wi
    dv_psixx = (
        -du_psixx_dx - dv_psixx_dy + 2.0 * omegalog * psixy + 2.0 * bxx + mxx / wi - tf_psixx
    )
    dv_psixy = (
        -du_psixy_dx
        - dv_psixy_dy
        + omegalog * (psiyy - psixx)
        + 2.0 * bxy
        + mxy / wi
        - tf_psixy
    )
    dv_psiyy = (
        -du_psiyy_dx - dv_psiyy_dy - 2.0 * omegalog * psixy + 2.0 * byy + myy / wi - tf_psiyy
    )

    return dvz, dv_psixx, dv_psixy, dv_psiyy


def euler_step_log_conformation(
    wz: FloatArray,
    psixx: FloatArray,
    psixy: FloatArray,
    psiyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfwz: FloatArray,
    tf_psixx: FloatArray,
    tf_psixy: FloatArray,
    tf_psiyy: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return ``(wz, Psixx, Psixy, Psiyy)`` advanced by one explicit Euler step.

    Port of the ``wz``/``Psixx``/``Psixy``/``Psiyy`` update inside
    ``subroutine solver_log_euler`` (``ty_sim == 22``, ``euler_rk.f``).
    Only the interior points are updated; the boundary is copied unchanged.
    """

    dvz, dv_psixx, dv_psixy, dv_psiyy = log_conformation_time_derivatives(
        wz, psixx, psixy, psiyy, ux, uy, tfwz, tf_psixx, tf_psixy, tf_psiyy, grid, parameters
    )
    dt = parameters.dt

    wz_new = np.array(wz, dtype=np.float64, copy=True)
    psixx_new = np.array(psixx, dtype=np.float64, copy=True)
    psixy_new = np.array(psixy, dtype=np.float64, copy=True)
    psiyy_new = np.array(psiyy, dtype=np.float64, copy=True)

    interior = (slice(1, -1), slice(1, -1))
    wz_new[interior] = wz[interior] + dt * dvz[interior]
    psixx_new[interior] = psixx[interior] + dt * dv_psixx[interior]
    psixy_new[interior] = psixy[interior] + dt * dv_psixy[interior]
    psiyy_new[interior] = psiyy[interior] + dt * dv_psiyy[interior]

    return wz_new, psixx_new, psixy_new, psiyy_new


__all__ = [
    "euler_step_log_conformation",
    "log_conformation_time_derivatives",
]
