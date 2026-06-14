"""One time-marching step for MMS Case 12 (Etapa 10).

This module combines the explicit Euler step from
:mod:`mms_viscoelastic.numerics.time_stepping` (Etapa 9) with the Poisson
solver from :mod:`mms_viscoelastic.numerics.poisson` (Etapa 8) to perform one
full time step of the Case 12 ``ty_sim == 1`` formulation:

1. Evaluate the manufactured source terms (``tfwz``, ``tfTxx``, ``tfTxy``,
   ``tfTyy``) at the current time ``t``.
2. Advance ``wz``, ``Txx``, ``Txy``, ``Tyy`` at the interior points by one
   explicit Euler step (:func:`euler_step_wz`, :func:`euler_step_stress`).
3. Re-impose Dirichlet boundary conditions at ``t + dt`` for ``wz``, ``Txx``,
   ``Txy``, ``Tyy``, ``psi``, ``u`` and ``v`` from the exact manufactured
   solution (:func:`evaluate_case_12_exact`).
4. Re-solve ``laplacian(psi) = wz`` for the interior of ``psi``
   (:func:`solve_poisson_dirichlet`), using the boundary values from step 3.
5. Recompute ``u = d(psi)/dy`` and ``v = -d(psi)/dx`` at the interior points
   from the updated ``psi``.

This is a verbatim-in-spirit port of the sequence
``call source_term(time + dt)`` / ``call boundaries_condictions(time + dt)``
/ ``call loop(2)`` in ``subroutine solver_euler``
(``src/fortran_original/4Order_LPOG_Euler/euler_rk.f``).

``wz`` boundary condition: a deliberate simplification
------------------------------------------------------
In the Fortran reference, ``boundaries_condictions_case_12`` re-imposes the
exact manufactured solution on the boundary for ``u``, ``v``, ``psi``,
``Txx``, ``Txy`` and ``Tyy``, but *not* for ``wz``. Instead it calls
``boundarie_vorticity_cavity_flow``
(``src/fortran_original/4Order_LPOG_Euler/boundarie_vorticity_cavity_flow.f``),
which derives the ``wz`` boundary values from one-sided 4th-order
finite-difference stencils applied to ``psi`` (bottom/left/right) and ``u``
(top) -- machinery tied to the 4th-order LPOG compact discretization.

Porting that stencil would reintroduce 4th-order-specific operators that
Etapas 8 and 9 deliberately left out of scope
(see ``docs/case-12-poisson-solver.md`` and
``docs/case-12-euler-step-comparison.md``). Instead, this module re-imposes
the ``wz`` boundary directly from :func:`evaluate_case_12_exact` at
``t + dt``, consistently with the other boundary fields. Both choices impose
the *exact* ``wz`` boundary value in the limit of grid refinement; the
Fortran one-sided stencil is simply a discrete approximation of it that is
consistent with the 4th-order interior operators. This is a deliberate
deviation from the Fortran reference and is documented in
``docs/case-12-time-marching-comparison.md``.

Multi-step marching (Etapa 11)
-------------------------------
:func:`case12_multistep_march` chains :func:`case12_time_step` over several
``dt`` increments, advancing the physical time by ``step_parameters.dt`` at
each iteration. See ``docs/case-12-multistep-marching-comparison.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.numerics.derivatives import derivative_x, derivative_y
from mms_viscoelastic.numerics.poisson import solve_poisson_dirichlet
from mms_viscoelastic.numerics.time_stepping import (
    EulerStepParameters,
    euler_step_stress,
    euler_step_wz,
)
from mms_viscoelastic.source_terms import Case12SourceTermParameters, evaluate_case_12_source_terms

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class TimeMarchingState:
    """Case 12 field state ``(u, v, wz, psi, Txx, Txy, Tyy)`` at a given time."""

    u: FloatArray
    v: FloatArray
    wz: FloatArray
    psi: FloatArray
    txx: FloatArray
    txy: FloatArray
    tyy: FloatArray


def _require_array(value: FloatArray | None, name: str) -> FloatArray:
    """Return ``value``, raising ``ValueError`` if it was not computed."""

    if value is None:
        raise ValueError(f"Field component {name!r} was not computed.")
    return value


def _apply_boundary(field: FloatArray, boundary_values: FloatArray) -> FloatArray:
    """Return ``field`` with its outer boundary overwritten by ``boundary_values``."""

    result = np.array(field, dtype=np.float64, copy=True)
    result[0, :] = boundary_values[0, :]
    result[-1, :] = boundary_values[-1, :]
    result[:, 0] = boundary_values[:, 0]
    result[:, -1] = boundary_values[:, -1]
    return result


def case12_time_step(
    state: TimeMarchingState,
    t: float,
    grid: Grid2D,
    step_parameters: EulerStepParameters | None = None,
    field_parameters: Case12Parameters | None = None,
    source_parameters: Case12SourceTermParameters | None = None,
) -> TimeMarchingState:
    """Advance ``state`` from time ``t`` to ``t + dt`` (one full Case 12 step).

    Parameters
    ----------
    state:
        Field state at time ``t``.
    t:
        Current physical time. The manufactured source terms are evaluated at
        ``t``, and the Dirichlet boundary conditions are re-imposed at
        ``t + step_parameters.dt``.
    grid:
        Uniform Cartesian grid.
    step_parameters, field_parameters, source_parameters:
        Defaults match :class:`EulerStepParameters`, :class:`Case12Parameters`
        and :class:`Case12SourceTermParameters`, which is also the parameter
        set used by Etapas 8 and 9.

    Returns
    -------
    TimeMarchingState
        Field state at time ``t + dt``.
    """
    step_params = step_parameters or EulerStepParameters()
    field_params = field_parameters or Case12Parameters()
    source_params = source_parameters or Case12SourceTermParameters()
    step_params.validate()

    source_terms = evaluate_case_12_source_terms(grid, t=t, parameters=source_params)
    tfwz = _require_array(source_terms.tfwz, "tfwz")
    tftxx = _require_array(source_terms.tfTxx, "tfTxx")
    tftxy = _require_array(source_terms.tfTxy, "tfTxy")
    tftyy = _require_array(source_terms.tfTyy, "tfTyy")

    wz_step = euler_step_wz(
        state.wz, state.txx, state.txy, state.tyy, state.u, state.v, tfwz, grid, step_params
    )
    txx_step, txy_step, tyy_step = euler_step_stress(
        state.txx, state.txy, state.tyy, state.u, state.v, tftxx, tftxy, tftyy, grid, step_params
    )

    t_next = t + step_params.dt
    exact_next = evaluate_case_12_exact(grid, t=t_next, parameters=field_params)
    wz_exact = _require_array(exact_next.wz, "wz")
    psi_exact = _require_array(exact_next.psi, "psi")
    txx_exact = _require_array(exact_next.txx, "txx")
    txy_exact = _require_array(exact_next.txy, "txy")
    tyy_exact = _require_array(exact_next.tyy, "tyy")
    u_exact = _require_array(exact_next.u, "u")
    v_exact = _require_array(exact_next.v, "v")

    wz_new = _apply_boundary(wz_step, wz_exact)
    txx_new = _apply_boundary(txx_step, txx_exact)
    txy_new = _apply_boundary(txy_step, txy_exact)
    tyy_new = _apply_boundary(tyy_step, tyy_exact)

    psi_new = solve_poisson_dirichlet(wz_new, grid, psi_exact)

    u_full = derivative_y(psi_new, grid.dy)
    v_full = -derivative_x(psi_new, grid.dx)
    u_new = _apply_boundary(u_full, u_exact)
    v_new = _apply_boundary(v_full, v_exact)

    return TimeMarchingState(
        u=u_new,
        v=v_new,
        wz=wz_new,
        psi=psi_new,
        txx=txx_new,
        txy=txy_new,
        tyy=tyy_new,
    )


def case12_multistep_march(
    state: TimeMarchingState,
    t0: float,
    n_steps: int,
    grid: Grid2D,
    step_parameters: EulerStepParameters | None = None,
    field_parameters: Case12Parameters | None = None,
    source_parameters: Case12SourceTermParameters | None = None,
) -> TimeMarchingState:
    """Advance ``state`` from ``t0`` by ``n_steps`` full time-marching steps.

    Repeatedly applies :func:`case12_time_step`, advancing the physical time
    by ``step_parameters.dt`` at each iteration.

    Parameters
    ----------
    state:
        Field state at time ``t0``.
    t0:
        Initial physical time.
    n_steps:
        Number of full time-marching steps to perform. Must be non-negative.
    grid:
        Uniform Cartesian grid.
    step_parameters, field_parameters, source_parameters:
        Forwarded to :func:`case12_time_step` at every step. Defaults match
        :class:`EulerStepParameters`, :class:`Case12Parameters` and
        :class:`Case12SourceTermParameters`.

    Returns
    -------
    TimeMarchingState
        Field state at time ``t0 + n_steps * step_parameters.dt``.
    """
    if n_steps < 0:
        raise ValueError(f"n_steps must be non-negative, got {n_steps}.")

    step_params = step_parameters or EulerStepParameters()
    step_params.validate()

    current = state
    for n in range(n_steps):
        current = case12_time_step(
            current,
            t0 + n * step_params.dt,
            grid,
            step_parameters=step_params,
            field_parameters=field_parameters,
            source_parameters=source_parameters,
        )
    return current


__all__ = [
    "TimeMarchingState",
    "case12_multistep_march",
    "case12_time_step",
]
