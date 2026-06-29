"""Explicit Euler time stepping for MMS Case 12 (``ty_sim == 1`` formulation).

This module ports the right-hand sides used by the Fortran subroutine
``drv`` (``src/fortran_original/4Order_LPOG_Euler/drvs.f``) and the explicit
Euler update of ``subroutine solver_euler``
(``src/fortran_original/4Order_LPOG_Euler/euler_rk.f``), for the governing
variables ``wz``, ``Txx``, ``Txy`` and ``Tyy``.

Sign convention
----------------
As documented in :mod:`mms_viscoelastic.numerics.poisson`, this package
follows ``case_12.f``'s convention ``wz = laplacian(psi)``. The vorticity
transport right-hand side below is the verbatim port of ``dvz`` from
``drv`` under that convention.

Discretization
--------------
As of Etapa 17, the derivatives are evaluated with the same 4th-order LPOG
compact (Padé-type) operators the Fortran ``drv`` actually uses
(``derx``/``dery``/``derxx``/``deryy`` in ``drvs.f`` are the compact
operators from ``derivs.f``, ported in Etapa 13 as
:func:`mms_viscoelastic.numerics.derivatives.compact_derivative_x`/``_y`` and
:func:`mms_viscoelastic.numerics.derivatives.compact_second_derivative_x`/``_y``).
Before Etapa 17, this module used the second-order operators as a
deliberate, documented scope reduction (Etapa 9); that scope reduction is
now closed, see ``docs/case-12-4th-order-time-marching.md``. The Lele
compact filter (``filter_trid_x``/``filter_trid_y`` in ``loop.f``) is a
separate step applied once per time step by
:mod:`mms_viscoelastic.numerics.time_marching`, matching ``subroutine loop``
in ``main.f`` -- it is not part of the right-hand sides in this module.

Only interior grid points (``[1:-1, 1:-1]``) are updated; boundary values are
copied unchanged from the input arrays. Re-imposing Dirichlet boundary
conditions, applying the Lele filter and recomputing ``psi``/``ux``/``uy``
from the updated ``wz`` (the Fortran ``boundaries_condictions`` and ``loop``
calls) is left to the multi-step time-marching solver.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.numerics.derivatives import (
    compact_derivative_x as derivative_x,
)
from mms_viscoelastic.numerics.derivatives import (
    compact_derivative_y as derivative_y,
)
from mms_viscoelastic.numerics.derivatives import (
    compact_second_derivative_x as second_derivative_x,
)
from mms_viscoelastic.numerics.derivatives import (
    compact_second_derivative_y as second_derivative_y,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class EulerStepParameters:
    """Physical and time-step parameters for the Case 12 Euler step.

    The defaults (``reynolds=1.0``, ``wi=1.0``, ``beta=0.5``) match
    :class:`mms_viscoelastic.manufactured_solutions.case_12.Case12Parameters`
    and :class:`mms_viscoelastic.source_terms.case_12.Case12SourceTermParameters`,
    which are the values already validated against the Fortran field and
    source-term references (FASE V-VII). ``dt=1.0e-5`` matches both
    ``src/fortran_original/4Order_LPOG_Euler/par.nn`` and the thesis figure
    filenames (``Dt_1e-05``).
    """

    reynolds: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    dt: float = 1.0e-5

    def validate(self) -> None:
        """Validate the parameter values used by the Euler-step formulas."""
        if self.reynolds == 0.0:
            raise ValueError("reynolds must be nonzero.")
        if self.wi == 0.0:
            raise ValueError("wi must be nonzero.")
        if self.beta == 1.0:
            raise ValueError("beta must be different from 1.0 (division by 1 - beta).")
        if self.dt <= 0.0:
            raise ValueError("dt must be positive.")


def wz_time_derivative(
    wz: FloatArray,
    txx: FloatArray,
    txy: FloatArray,
    tyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfwz: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> FloatArray:
    """Return ``dvz``, the right-hand side of the ``wz`` transport equation.

    Port of ``dvz`` in ``drv`` (``drvs.f``):

    ``dvz = -d(ux*wz)/dx - d(uy*wz)/dy + (beta/Re) * laplacian(wz)``
    ``     + d2Txx/dxdy + d2Txy/dy2 - d2Txy/dx2 - d2Tyy/dxdy - tfwz``
    """
    parameters.validate()
    dx, dy = grid.dx, grid.dy

    duwzdx = derivative_x(ux * wz, dx)
    dvwzdy = derivative_y(uy * wz, dy)
    diffusion = second_derivative_x(wz, dx) + second_derivative_y(wz, dy)

    d2txx_dxdy = derivative_y(derivative_x(txx, dx), dy)
    d2txy_dx2 = second_derivative_x(txy, dx)
    d2txy_dy2 = second_derivative_y(txy, dy)
    d2tyy_dxdy = derivative_y(derivative_x(tyy, dx), dy)

    return (
        -duwzdx
        - dvwzdy
        + (parameters.beta / parameters.reynolds) * diffusion
        + d2txx_dxdy
        + d2txy_dy2
        - d2txy_dx2
        - d2tyy_dxdy
        - tfwz
    )


def stress_time_derivatives(
    txx: FloatArray,
    txy: FloatArray,
    tyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfTxx: FloatArray,
    tfTxy: FloatArray,
    tfTyy: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return ``(dvTxx, dvTxy, dvTyy)``, the UCM/Oldroyd-B/Giesekus right-hand sides.

    Port of ``dvTxx``, ``dvTxy`` and ``dvTyy`` in ``drv`` (``drvs.f``).
    """
    parameters.validate()
    p = parameters
    dx, dy = grid.dx, grid.dy

    dudx = derivative_x(ux, dx)
    dudy = derivative_y(ux, dy)
    dvdx = derivative_x(uy, dx)
    dvdy = derivative_y(uy, dy)

    duTxxdx = derivative_x(ux * txx, dx)
    dvTxxdy = derivative_y(uy * txx, dy)
    duTxydx = derivative_x(ux * txy, dx)
    dvTxydy = derivative_y(uy * txy, dy)
    duTyydx = derivative_x(ux * tyy, dx)
    dvTyydy = derivative_y(uy * tyy, dy)

    one_minus_beta = 1.0 - p.beta
    f_tij = 1.0 + (p.epsilon * p.reynolds * p.wi / one_minus_beta) * (txx + tyy)

    dv_txx = (
        -f_tij * txx / p.wi
        - duTxxdx
        - dvTxxdy
        + 2.0 * txx * dudx
        + 2.0 * txy * dudy
        - 2.0 * p.xi * txx * dudx
        - p.xi * txy * (dudy + dvdx)
        - (p.alpha_g * p.reynolds / one_minus_beta) * (txx * txx + txy * txy)
        + (2.0 * one_minus_beta / (p.reynolds * p.wi)) * dudx
        + tfTxx / p.wi
    )

    dv_txy = (
        -f_tij * txy / p.wi
        - duTxydx
        - dvTxydy
        + txx * dvdx
        + tyy * dudy
        - p.xi * txy * dudx
        - 0.5 * p.xi * tyy * (dudy + dvdx)
        - 0.5 * p.xi * txx * (dudy + dvdx)
        - p.xi * txy * dvdy
        - (p.alpha_g * p.reynolds / one_minus_beta) * (txx * txy + txy * tyy)
        + (one_minus_beta / (p.reynolds * p.wi)) * (dudy + dvdx)
        + tfTxy / p.wi
    )

    dv_tyy = (
        -f_tij * tyy / p.wi
        - duTyydx
        - dvTyydy
        + 2.0 * txy * dvdx
        + 2.0 * tyy * dvdy
        - p.xi * txy * (dvdx + dudy)
        - 2.0 * p.xi * tyy * dvdy
        - (p.alpha_g * p.reynolds / one_minus_beta) * (txy * txy + tyy * tyy)
        + (2.0 * one_minus_beta / (p.reynolds * p.wi)) * dvdy
        + tfTyy / p.wi
    )

    return dv_txx, dv_txy, dv_tyy


def _euler_update_interior(field: FloatArray, time_derivative: FloatArray, dt: float) -> FloatArray:
    """Return ``field`` with interior points advanced by ``time_derivative * dt``."""
    updated = np.array(field, dtype=np.float64, copy=True)
    updated[1:-1, 1:-1] = field[1:-1, 1:-1] + time_derivative[1:-1, 1:-1] * dt
    return updated


def euler_step_wz(
    wz: FloatArray,
    txx: FloatArray,
    txy: FloatArray,
    tyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfwz: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> FloatArray:
    """Return ``wz`` advanced by one explicit Euler step (interior points only)."""
    dvz = wz_time_derivative(wz, txx, txy, tyy, ux, uy, tfwz, grid, parameters)
    return _euler_update_interior(wz, dvz, parameters.dt)


def euler_step_stress(
    txx: FloatArray,
    txy: FloatArray,
    tyy: FloatArray,
    ux: FloatArray,
    uy: FloatArray,
    tfTxx: FloatArray,
    tfTxy: FloatArray,
    tfTyy: FloatArray,
    grid: Grid2D,
    parameters: EulerStepParameters,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return ``(Txx, Txy, Tyy)`` advanced by one explicit Euler step (interior points only)."""
    dv_txx, dv_txy, dv_tyy = stress_time_derivatives(
        txx, txy, tyy, ux, uy, tfTxx, tfTxy, tfTyy, grid, parameters
    )
    return (
        _euler_update_interior(txx, dv_txx, parameters.dt),
        _euler_update_interior(txy, dv_txy, parameters.dt),
        _euler_update_interior(tyy, dv_tyy, parameters.dt),
    )


__all__ = [
    "EulerStepParameters",
    "euler_step_stress",
    "euler_step_wz",
    "stress_time_derivatives",
    "wz_time_derivative",
]
