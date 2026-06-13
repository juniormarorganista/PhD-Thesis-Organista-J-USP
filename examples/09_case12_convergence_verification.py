"""Grid-convergence verification of the finite-difference operators (case 12).

This script applies the second-order operators ``derivative_x``,
``derivative_y`` to the exact manufactured fields of case 12 and checks that
the discrete stream-function relations

    u = d(psi)/dy
    v = -d(psi)/dx
    wz = d(u)/dy - d(v)/dx

are reproduced with the expected order of accuracy ``p = 2`` as the grid is
refined.

Run from the repository root after installing the package in editable mode:

    python examples/09_case12_convergence_verification.py

Optional custom grid sizes (number of cells per direction):

    python examples/09_case12_convergence_verification.py --grid-sizes 8 16 32 64 128
"""

from __future__ import annotations

import argparse

from mms_viscoelastic.numerics import QUANTITIES, run_convergence_study


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Grid-convergence study for MMS case 12 finite-difference operators."
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
    report = run_convergence_study(args.grid_sizes, t=args.time)

    print("MMS case 12 - grid-convergence verification")
    print("Relations checked: u = d(psi)/dy, v = -d(psi)/dx, wz = d(u)/dy - d(v)/dx\n")

    header = f"{'n':>5s} {'h':>12s}"
    for quantity in QUANTITIES:
        header += f" {('L2_' + quantity):>14s} {('Linf_' + quantity):>14s}"
    print(header)

    for level in report.levels:
        row = f"{level.n:5d} {level.h:12.6e}"
        for quantity in QUANTITIES:
            row += f" {level.l2_errors[quantity]:14.6e} {level.linf_errors[quantity]:14.6e}"
        print(row)

    print("\nEstimated order of accuracy (log-log fit across all grids):")
    for quantity in QUANTITIES:
        print(
            f"  {quantity:>3s}: L2 order = {report.l2_orders[quantity]:.4f}, "
            f"Linf order = {report.linf_orders[quantity]:.4f}"
        )


if __name__ == "__main__":
    main()
