"""Direct solver for the 2D Poisson equation with Dirichlet boundary conditions.

This module solves ``laplacian(phi) = source`` on a uniform Cartesian grid
using a sparse direct solve over the standard second-order, five-point
finite-difference Laplacian. It is the missing piece needed to recompute the
stream function ``psi`` from the vorticity ``wz`` (``laplacian(psi) = wz`` in
this package's sign convention, see :mod:`mms_viscoelastic.numerics.convergence`)
during time-marching.

The original Fortran solver (``src/fortran_original/4Order_LPOG_Euler/poisson.f``)
uses a multigrid V-cycle with SOR smoothing. For the manufactured-solution
verification grids used here (up to roughly 65 x 65 points), a sparse direct
solve is equivalent in accuracy and much simpler to validate; the multigrid
algorithm itself is not ported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, diags, eye, kron
from scipy.sparse.linalg import spsolve

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions.case_12 import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.numerics.errors import estimate_order, l2_error, linf_error

FloatArray = NDArray[np.float64]


def _second_difference_matrix(n: int, h2: float) -> csr_matrix:
    """Return the ``n x n`` 1D second-difference matrix divided by ``h2``."""
    if n == 1:
        return csr_matrix(np.array([[-2.0 / h2]]))
    return diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)).tocsr() / h2


def solve_poisson_dirichlet(
    source: FloatArray, grid: Grid2D, boundary_values: FloatArray
) -> FloatArray:
    """Solve ``laplacian(phi) = source`` with Dirichlet boundary conditions.

    Parameters
    ----------
    source:
        Right-hand side ``source[i, j]``, shape ``(grid.ny, grid.nx)``. Only
        the interior entries (excluding the domain boundary) are used.
    grid:
        Uniform Cartesian grid.
    boundary_values:
        Field values to impose on the domain boundary, shape
        ``(grid.ny, grid.nx)``. Only the entries on the first/last row and
        column are used.

    Returns
    -------
    FloatArray
        Solution ``phi`` with shape ``(grid.ny, grid.nx)``: the boundary
        entries are copied from ``boundary_values`` and the interior entries
        solve the discrete Poisson equation exactly (to direct-solve
        round-off).

    Raises
    ------
    ValueError
        If ``source`` or ``boundary_values`` do not have shape
        ``(grid.ny, grid.nx)``, or if the grid has no interior points.
    """
    source = np.asarray(source, dtype=np.float64)
    boundary_values = np.asarray(boundary_values, dtype=np.float64)
    expected_shape = (grid.ny, grid.nx)
    if source.shape != expected_shape:
        raise ValueError(f"source must have shape {expected_shape}, got {source.shape}.")
    if boundary_values.shape != expected_shape:
        raise ValueError(
            f"boundary_values must have shape {expected_shape}, got {boundary_values.shape}."
        )

    n_interior_y = grid.ny - 2
    n_interior_x = grid.nx - 2
    if n_interior_y < 1 or n_interior_x < 1:
        raise ValueError("Grid must have at least one interior point in each direction.")

    dx2 = grid.dx**2
    dy2 = grid.dy**2

    dxx = _second_difference_matrix(n_interior_x, dx2)
    dyy = _second_difference_matrix(n_interior_y, dy2)
    laplacian_matrix = kron(eye(n_interior_y), dxx) + kron(dyy, eye(n_interior_x))

    boundary_term = np.zeros((n_interior_y, n_interior_x))
    boundary_term[:, 0] += boundary_values[1:-1, 0] / dx2
    boundary_term[:, -1] += boundary_values[1:-1, -1] / dx2
    boundary_term[0, :] += boundary_values[0, 1:-1] / dy2
    boundary_term[-1, :] += boundary_values[-1, 1:-1] / dy2

    rhs = source[1:-1, 1:-1].ravel() - boundary_term.ravel()
    phi_interior = spsolve(csr_matrix(laplacian_matrix), rhs)

    phi = boundary_values.copy()
    phi[1:-1, 1:-1] = phi_interior.reshape(n_interior_y, n_interior_x)
    return phi


@dataclass(frozen=True)
class PoissonConvergenceLevel:
    """Errors of a single grid in :func:`run_case12_poisson_convergence_study`."""

    n: int
    h: float
    l2_error: float
    linf_error: float


@dataclass(frozen=True)
class PoissonConvergenceReport:
    """Result of :func:`run_case12_poisson_convergence_study`."""

    levels: list[PoissonConvergenceLevel]
    l2_order: float
    linf_order: float

    def to_rows(self) -> list[dict[str, float | int]]:
        """Flatten the report into a list of dictionaries (one per grid)."""
        return [
            {
                "n": level.n,
                "h": level.h,
                "l2_error": level.l2_error,
                "linf_error": level.linf_error,
            }
            for level in self.levels
        ]


def run_case12_poisson_convergence_study(
    grid_sizes: list[int],
    *,
    t: float = 0.0,
    parameters: Case12Parameters | None = None,
) -> PoissonConvergenceReport:
    """Verify :func:`solve_poisson_dirichlet` against the case 12 MMS.

    For each grid, ``wz_exact`` is used as the Poisson source and
    ``psi_exact`` supplies the Dirichlet boundary values; the resulting
    ``psi_numerical`` is compared against ``psi_exact`` over the whole grid.

    Parameters
    ----------
    grid_sizes:
        Numbers of cells ``n`` along each direction. Each grid has
        ``nx = ny = n + 1`` points and mesh size ``h = 1 / n``. Must contain
        at least two values.
    t:
        Physical time at which the exact fields are evaluated.
    parameters:
        Case-12 parameters. Defaults to :class:`Case12Parameters` defaults.

    Returns
    -------
    PoissonConvergenceReport
        Per-grid errors and the estimated order of accuracy.

    Raises
    ------
    ValueError
        If ``grid_sizes`` has fewer than two elements.
    """
    if len(grid_sizes) < 2:
        raise ValueError("grid_sizes must contain at least two grid resolutions.")

    levels: list[PoissonConvergenceLevel] = []
    for n in grid_sizes:
        grid = Grid2D.from_bounds(nx=n + 1, ny=n + 1)
        field = evaluate_case_12_exact(grid, t=t, parameters=parameters)
        assert field.wz is not None
        assert field.psi is not None

        psi_numerical = solve_poisson_dirichlet(field.wz, grid, field.psi)
        levels.append(
            PoissonConvergenceLevel(
                n=n,
                h=grid.dx,
                l2_error=l2_error(psi_numerical, field.psi),
                linf_error=linf_error(psi_numerical, field.psi),
            )
        )

    h_values = np.array([level.h for level in levels], dtype=np.float64)
    l2_order = estimate_order(h_values, np.array([level.l2_error for level in levels]))
    linf_order = estimate_order(h_values, np.array([level.linf_error for level in levels]))

    return PoissonConvergenceReport(levels=levels, l2_order=l2_order, linf_order=linf_order)


__all__ = [
    "PoissonConvergenceLevel",
    "PoissonConvergenceReport",
    "run_case12_poisson_convergence_study",
    "solve_poisson_dirichlet",
]
