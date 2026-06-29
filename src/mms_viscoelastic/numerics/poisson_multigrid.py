"""Multigrid V-cycle Poisson solver, 4th-order compact (Mehrstellen) stencil.

Ports ``subroutine poi_psi4a`` and its dependents (``verifica_poi``,
``function_source``, ``sor``, ``defect``, ``restrict``, ``intpol``,
``correction``, ``tridv``) in
``src/fortran_original/4Order_LPOG_Euler/poisson.f`` (Etapa 17).

This is the Poisson solver actually used by the original time-marching loop
(``call poi_psi4a`` in ``subroutine loop``, ``main.f``). The pentadiagonal
``band5_poi``/``banbk5_poi`` routines in ``solvers.f`` are **not** this
solver: a repository-wide search shows they are never called anywhere in
the original Fortran source. They are dead code and are not ported.

Algorithm
---------
``poi_psi4a`` solves ``laplacian(psi) = wz`` (this package's sign
convention, see :mod:`mms_viscoelastic.numerics.convergence`) with a 4th-order
compact 9-point (Mehrstellen-type) discretization:

    (1+r)(NE+NW+SE+SW) + (10-2r)(N+S) + (-2+10r)(E+W) - 20(1+r)*P = dy2 * B(wz)

where ``r = dy2/dx2`` and ``B(wz) = wz_N + wz_S + wz_E + wz_W + 8*wz_P`` is a
compact 9-point weighting of the source term (``function_source``). This
linear system is solved by a **FAS (Full Approximation Scheme) multigrid
V-cycle**: line-relaxation SOR per mesh level (``sor``, solving a
tridiagonal system per grid line via the Thomas algorithm, ``tridv`` /
:func:`mms_viscoelastic.numerics.derivatives.solve_tridiagonal`), straight
injection + full-weighting restriction of the defect to a coarser mesh
(``restrict``), and bilinear-interpolation prolongation of the coarse-grid
correction back to the finer mesh (``intpol``, ``correction``).

Each mesh level halves the grid via ``n -> (n + 1) // 2``, matching
``par.nn``'s ``imax = jmax = 33``, ``msh = 4`` (levels of 33, 17, 9, 5
points). The line-relaxation sweep is **not** vectorizable across grid
lines: each line's right-hand side reads the *already-updated* neighboring
line (one direction) and the *not-yet-updated* one (the other direction),
i.e. it is a Gauss-Seidel-type zebra sweep, not Jacobi. This module loops
over grid lines in Python to reproduce that dependency exactly, matching
the Fortran reference to machine precision (see
``validation/fortran_reference/case12_poisson_multigrid/``).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.numerics.derivatives import solve_tridiagonal

FloatArray = NDArray[np.float64]

DEFAULT_N_LEVELS = 4
DEFAULT_MAX_VCYCLES = 81
DEFAULT_ERROR_TOL = 1.0e-6
_RELAXATION_FACTOR = 1.2
_RELAXATION_FACTOR_SINGLE_SWEEP = 1.0


def _function_source(
    wz: FloatArray, psi: FloatArray, n_levels: int, nx: int, ny: int
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Port of ``function_source``. Arrays are ``(level, i, j)`` (Fortran order)."""

    fc = np.zeros((n_levels, nx, ny), dtype=np.float64)
    sc = np.zeros((n_levels, nx, ny), dtype=np.float64)
    fct = np.zeros((n_levels, nx, ny), dtype=np.float64)

    fc[0] = psi
    sc[0, 1:-1, 1:-1] = (
        wz[1:-1, 0:-2] + wz[1:-1, 2:] + wz[2:, 1:-1] + wz[0:-2, 1:-1] + 8.0 * wz[1:-1, 1:-1]
    )
    return fc, sc, fct


def _sor(
    fc: FloatArray,
    sc: FloatArray,
    nptsx: int,
    nptsy: int,
    dy2: float,
    dyypdxx: float,
    itera: int,
    mesh: int,
) -> None:
    """Port of ``sor``: in-place line-relaxation smoothing of ``fc[mesh]``."""

    a = np.zeros(nptsy, dtype=np.float64)
    b = np.zeros(nptsy, dtype=np.float64)
    c = np.zeros(nptsy, dtype=np.float64)
    b[0] = 1.0
    b[-1] = 1.0
    a[1:-1] = 10.0 - 2.0 * dyypdxx
    b[1:-1] = -20.0 * (1.0 + dyypdxx)
    c[1:-1] = 10.0 - 2.0 * dyypdxx

    rf = _RELAXATION_FACTOR_SINGLE_SWEEP if itera == 1 else _RELAXATION_FACTOR
    fcm = fc[mesh, :nptsx, :nptsy]
    scm = sc[mesh, :nptsx, :nptsy]
    rhs = np.empty(nptsy, dtype=np.float64)

    for _ in range(itera):
        for i in range(1, nptsx - 1):
            rhs[0] = fcm[i, 0]
            rhs[1:-1] = scm[i, 1:-1] * dy2 - (
                (1.0 + dyypdxx)
                * (fcm[i + 1, 2:] + fcm[i - 1, 2:] + fcm[i + 1, 0:-2] + fcm[i - 1, 0:-2])
                + (-2.0 + 10.0 * dyypdxx) * (fcm[i + 1, 1:-1] + fcm[i - 1, 1:-1])
            )
            rhs[-1] = fcm[i, -1]
            solved = solve_tridiagonal(a, b, c, rhs[:, None])[:, 0]
            fcm[i, 1:-1] = fcm[i, 1:-1] + rf * (solved[1:-1] - fcm[i, 1:-1])


def _defect(
    fc: FloatArray,
    sc: FloatArray,
    fct: FloatArray,
    nptsx: int,
    nptsy: int,
    dx2: float,
    dy2: float,
    si: float,
    mesh: int,
) -> None:
    """Port of ``defect``: residual (``si < 0``) or FAS coarse RHS (``si > 0``)."""

    fcm = fc[mesh, :nptsx, :nptsy]
    b = (
        fcm[2:, 2:]
        + fcm[2:, 0:-2]
        + fcm[0:-2, 2:]
        + fcm[0:-2, 0:-2]
        + 10.0 * (fcm[2:, 1:-1] + fcm[0:-2, 1:-1])
        - 2.0 * (fcm[1:-1, 2:] + fcm[1:-1, 0:-2])
        - 20.0 * fcm[1:-1, 1:-1]
    )
    c = (
        fcm[2:, 2:]
        + fcm[2:, 0:-2]
        + fcm[0:-2, 2:]
        + fcm[0:-2, 0:-2]
        + 10.0 * (fcm[1:-1, 2:] + fcm[1:-1, 0:-2])
        - 2.0 * (fcm[2:, 1:-1] + fcm[0:-2, 1:-1])
        - 20.0 * fcm[1:-1, 1:-1]
    )
    fct[0, 1 : nptsx - 1, 1 : nptsy - 1] = sc[mesh, 1 : nptsx - 1, 1 : nptsy - 1] + si * (
        b / dx2 + c / dy2
    )

    if si < 0:
        return
    sc[mesh, 1 : nptsx - 1, 1 : nptsy - 1] = fct[0, 1 : nptsx - 1, 1 : nptsy - 1]


def _restrict(
    fc: FloatArray,
    sc: FloatArray,
    fct: FloatArray,
    nptsx: int,
    nptsy: int,
    dx2: float,
    dy2: float,
    mesh: int,
) -> tuple[int, int, float, float]:
    """Port of ``restrict``. Returns the updated ``(nptsx, nptsy, dx2, dy2)``."""

    new_nptsx = (nptsx + 1) // 2
    new_nptsy = (nptsy + 1) // 2
    new_dx2 = 4.0 * dx2
    new_dy2 = 4.0 * dy2

    fine = 2 + 2 * np.arange(new_nptsx - 2)  # fi = 3,5,... (0-based: 2,4,...)
    fine_j = 2 + 2 * np.arange(new_nptsy - 2)
    fi_grid, fj_grid = np.meshgrid(fine, fine_j, indexing="ij")

    injected = fc[mesh][fi_grid, fj_grid]
    fc[mesh + 1, 1 : new_nptsx - 1, 1 : new_nptsy - 1] = injected
    fct[mesh + 1, 1 : new_nptsx - 1, 1 : new_nptsy - 1] = injected

    fct0 = fct[0]
    sc[mesh + 1, 1 : new_nptsx - 1, 1 : new_nptsy - 1] = 6.25e-2 * (
        2.0
        * (
            fct0[fi_grid - 1, fj_grid]
            + fct0[fi_grid + 1, fj_grid]
            + fct0[fi_grid, fj_grid - 1]
            + fct0[fi_grid, fj_grid + 1]
        )
        + 4.0 * fct0[fi_grid, fj_grid]
        + fct0[fi_grid + 1, fj_grid + 1]
        + fct0[fi_grid + 1, fj_grid - 1]
        + fct0[fi_grid - 1, fj_grid + 1]
        + fct0[fi_grid - 1, fj_grid - 1]
    )

    return new_nptsx, new_nptsy, new_dx2, new_dy2


def _intpol(
    fc: FloatArray, nptsx: int, nptsy: int, dx2: float, dy2: float, mesh: int
) -> tuple[int, int, float, float]:
    """Port of ``intpol``. Returns the updated ``(nptsx, nptsy, dx2, dy2)``."""

    nx_total, ny_total = fc.shape[1], fc.shape[2]
    new_dx2 = 0.25 * dx2
    new_dy2 = 0.25 * dy2

    temp = np.zeros((nx_total, ny_total), dtype=np.float64)
    temp[0 : 2 * nptsx - 1 : 2, 0 : 2 * nptsy - 1 : 2] = fc[mesh, 0:nptsx, 0:nptsy]

    new_nptsx = 2 * nptsx - 1
    new_nptsy = 2 * nptsy - 1

    temp[1 : new_nptsx - 1 : 2, 0:new_nptsy:2] = 0.5 * (
        temp[0 : new_nptsx - 2 : 2, 0:new_nptsy:2] + temp[2:new_nptsx:2, 0:new_nptsy:2]
    )
    temp[0:new_nptsx, 1 : new_nptsy - 1 : 2] = 0.5 * (
        temp[0:new_nptsx, 2:new_nptsy:2] + temp[0:new_nptsx, 0 : new_nptsy - 2 : 2]
    )

    fc[mesh - 1, 1 : new_nptsx - 1, 1 : new_nptsy - 1] += temp[
        1 : new_nptsx - 1, 1 : new_nptsy - 1
    ]

    return new_nptsx, new_nptsy, new_dx2, new_dy2


def _correction(fc: FloatArray, fct: FloatArray, nptsx: int, nptsy: int, mesh: int) -> None:
    """Port of ``correction``."""

    fc[mesh, 1 : nptsx - 1, 1 : nptsy - 1] -= fct[mesh, 1 : nptsx - 1, 1 : nptsy - 1]


def solve_poisson_multigrid_4th_order(
    wz: FloatArray,
    grid: Grid2D,
    psi_boundary: FloatArray,
    *,
    n_levels: int = DEFAULT_N_LEVELS,
    max_vcycles: int = DEFAULT_MAX_VCYCLES,
    error_tol: float = DEFAULT_ERROR_TOL,
) -> FloatArray:
    """Solve ``laplacian(psi) = wz`` with the 4th-order compact multigrid V-cycle.

    Port of ``poi_psi4a`` (``src/fortran_original/4Order_LPOG_Euler/poisson.f``).

    Parameters
    ----------
    wz:
        Right-hand side, shape ``(grid.ny, grid.nx)``.
    grid:
        Uniform Cartesian grid.
    psi_boundary:
        Field used both as the Dirichlet boundary condition (only the outer
        row/column are used, matching the Fortran: ``sor``'s tridiagonal
        system pins row 1 and row ``nptsy`` to whatever ``fc`` currently
        holds there, which is never reassigned) and as the V-cycle's initial
        guess for the interior (``function_source``'s ``fc(1,:,:) = psi``).
        Shape ``(grid.ny, grid.nx)``.
    n_levels:
        Number of multigrid levels. Defaults to ``4``, matching ``par.nn``'s
        ``msh = 4``. ``grid.nx``/``grid.ny`` must support ``n_levels - 1``
        halvings of the form ``n -> (n + 1) // 2`` while staying ``>= 3``.
    max_vcycles:
        Maximum number of V-cycles. Defaults to ``81``, matching the
        Fortran ``do vcycles = 1, 81`` loop bound.
    error_tol:
        Convergence tolerance on the max-abs finest-grid residual.
        Defaults to ``1e-6``, matching ``erromax`` in ``poi_psi4a``.

    Returns
    -------
    FloatArray
        ``psi`` with shape ``(grid.ny, grid.nx)``.
    """

    wz = np.asarray(wz, dtype=np.float64)
    psi_boundary = np.asarray(psi_boundary, dtype=np.float64)
    expected_shape = (grid.ny, grid.nx)
    if wz.shape != expected_shape:
        raise ValueError(f"wz must have shape {expected_shape}, got {wz.shape}.")
    if psi_boundary.shape != expected_shape:
        raise ValueError(
            f"psi_boundary must have shape {expected_shape}, got {psi_boundary.shape}."
        )
    if n_levels < 1:
        raise ValueError("n_levels must be at least 1.")

    # Transpose to Fortran (i, j) = (x, y) array order for literal index fidelity.
    wz_ij = wz.T
    psi_ij = psi_boundary.T
    nx, ny = wz_ij.shape

    fc, sc, fct = _function_source(wz_ij, psi_ij, n_levels, nx, ny)

    nptsx, nptsy = nx, ny
    dx2, dy2 = grid.dx**2, grid.dy**2
    itera = 2

    for _ in range(max_vcycles):
        for mesh in range(n_levels - 1):
            dyypdxx = dy2 / dx2
            _sor(fc, sc, nptsx, nptsy, dy2, dyypdxx, itera, mesh)
            _defect(fc, sc, fct, nptsx, nptsy, dx2, dy2, -1.0, mesh)
            if mesh == 0:
                err = np.max(np.abs(fct[0, 1 : nptsx - 1, 1 : nptsy - 1]))
                if err < error_tol:
                    return fc[0].T.copy()
            nptsx, nptsy, dx2, dy2 = _restrict(fc, sc, fct, nptsx, nptsy, dx2, dy2, mesh)
            dyypdxx = dy2 / dx2
            _defect(fc, sc, fct, nptsx, nptsy, dx2, dy2, 1.0, mesh + 1)

        dyypdxx = dy2 / dx2
        _sor(fc, sc, nptsx, nptsy, dy2, dyypdxx, 100 * itera, n_levels - 1)

        for mesh in range(n_levels - 1, 0, -1):
            _correction(fc, fct, nptsx, nptsy, mesh)
            nptsx, nptsy, dx2, dy2 = _intpol(fc, nptsx, nptsy, dx2, dy2, mesh)
            dyypdxx = dy2 / dx2
            _sor(fc, sc, nptsx, nptsy, dy2, dyypdxx, 1, mesh - 1)

    return fc[0].T.copy()


__all__ = [
    "DEFAULT_ERROR_TOL",
    "DEFAULT_MAX_VCYCLES",
    "DEFAULT_N_LEVELS",
    "solve_poisson_multigrid_4th_order",
]
