"""Grid-convergence verification of the direct Poisson solver (case 12).

This script solves ``laplacian(psi) = wz`` with Dirichlet boundary conditions
taken from the exact ``psi`` of MMS case 12 (``solve_poisson_dirichlet``), and
checks that the recovered ``psi`` converges to the exact field at second order
as the grid is refined.

Run from the repository root after installing the package in editable mode:

    python examples/12_case12_poisson_solver.py

Optional custom grid sizes (number of cells per direction):

    python examples/12_case12_poisson_solver.py --grid-sizes 8 16 32 64 128
"""

from __future__ import annotations

import argparse

from mms_viscoelastic.numerics import run_case12_poisson_convergence_study


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Grid-convergence study for the case 12 Poisson solver (laplacian(psi) = wz)."
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
    report = run_case12_poisson_convergence_study(args.grid_sizes, t=args.time)

    print("MMS case 12 - Poisson solver verification")
    print("Relation checked: laplacian(psi) = wz, Dirichlet BC = psi_exact\n")

    print(f"{'n':>5s} {'h':>12s} {'L2_psi':>14s} {'Linf_psi':>14s}")
    for level in report.levels:
        print(f"{level.n:5d} {level.h:12.6e} {level.l2_error:14.6e} {level.linf_error:14.6e}")

    print("\nEstimated order of accuracy (log-log fit across all grids):")
    print(f"  psi: L2 order = {report.l2_order:.4f}, Linf order = {report.linf_order:.4f}")


if __name__ == "__main__":
    main()
