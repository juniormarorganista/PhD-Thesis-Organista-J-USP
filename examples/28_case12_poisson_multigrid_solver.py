"""Grid-convergence verification of the 4th-order multigrid Poisson solver (Etapa 17).

This script solves ``laplacian(psi) = wz`` with the FAS multigrid V-cycle
(``solve_poisson_multigrid_4th_order``, a port of ``poi_psi4a`` in
``src/fortran_original/4Order_LPOG_Euler/poisson.f``), using Dirichlet
boundary conditions taken from the exact ``psi`` of MMS case 12, and checks
that the recovered ``psi`` converges to the exact field at fourth order as
the grid is refined -- unlike the direct second-order solver in
``examples/12_case12_poisson_solver.py``.

Run from the repository root after installing the package in editable mode:

    python examples/28_case12_poisson_multigrid_solver.py

Optional custom grid sizes (number of cells per direction):

    python examples/28_case12_poisson_multigrid_solver.py --grid-sizes 8 16 32 64
"""

from __future__ import annotations

import argparse

import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.numerics.errors import estimate_order, l2_error, linf_error
from mms_viscoelastic.numerics.poisson_multigrid import solve_poisson_multigrid_4th_order


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Grid-convergence study for the case 12 multigrid Poisson solver."
    )
    parser.add_argument(
        "--grid-sizes",
        type=int,
        nargs="+",
        default=[8, 16, 32, 64],
        help="Number of cells per direction for each refinement level.",
    )
    parser.add_argument("--time", type=float, default=0.0, help="Physical time.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("MMS case 12 - 4th-order compact multigrid Poisson solver verification")
    print("Relation checked: laplacian(psi) = wz, Dirichlet BC = psi_exact\n")
    print(f"{'n':>5s} {'h':>12s} {'L2_psi':>14s} {'Linf_psi':>14s}")

    l2_errors = []
    linf_errors = []
    h_values = []
    for n in args.grid_sizes:
        grid = Grid2D.from_bounds(nx=n + 1, ny=n + 1)
        fields = evaluate_case_12_exact(grid, t=args.time, parameters=Case12Parameters())
        assert fields.wz is not None
        assert fields.psi is not None

        psi_numerical = solve_poisson_multigrid_4th_order(fields.wz, grid, fields.psi)
        l2 = l2_error(psi_numerical, fields.psi)
        linf = linf_error(psi_numerical, fields.psi)
        l2_errors.append(l2)
        linf_errors.append(linf)
        h_values.append(grid.dx)
        print(f"{n:5d} {grid.dx:12.6e} {l2:14.6e} {linf:14.6e}")

    l2_order = estimate_order(np.array(h_values), np.array(l2_errors))
    linf_order = estimate_order(np.array(h_values), np.array(linf_errors))
    print("\nEstimated order of accuracy (log-log fit across all grids):")
    print(f"  psi: L2 order = {l2_order:.4f}, Linf order = {linf_order:.4f}")


if __name__ == "__main__":
    main()
