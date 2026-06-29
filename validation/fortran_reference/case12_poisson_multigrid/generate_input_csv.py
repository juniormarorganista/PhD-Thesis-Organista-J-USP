"""Generate the input CSV consumed by the multigrid Poisson Fortran reference.

Evaluates the already-validated Python Case 12 exact fields
(``evaluate_case_12_exact``) at a given time on a uniform grid and writes::

    x,y,wz,psi

The Fortran reference program (``case12_poisson_multigrid_reference.f90``)
reads this CSV, solves ``laplacian(psi) = wz`` via the FAS multigrid V-cycle
(``poi_psi4a``), using the input ``psi`` both as the Dirichlet boundary
condition and as the initial guess for the V-cycle (matching
``function_source``'s ``fc(1,:,:) = psi``), and writes the post-solve
``psi``.

Usage::

    python validation/fortran_reference/case12_poisson_multigrid/generate_input_csv.py \\
        validation/reference_data/case12_poisson_multigrid_input_33x33.csv 33 33
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact

DEFAULT_OUTPUT = "validation/reference_data/case12_poisson_multigrid_input_33x33.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", default=DEFAULT_OUTPUT)
    parser.add_argument("nx", nargs="?", type=int, default=33)
    parser.add_argument("ny", nargs="?", type=int, default=33)
    parser.add_argument("t", nargs="?", type=float, default=0.0)
    args = parser.parse_args()

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    fields = evaluate_case_12_exact(grid, t=args.t, parameters=Case12Parameters())

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "y", "wz", "psi"])
        for j in range(grid.ny):
            for i in range(grid.nx):
                writer.writerow(
                    [
                        f"{grid.X[j, i]:.17e}",
                        f"{grid.Y[j, i]:.17e}",
                        f"{fields.wz[j, i]:.17e}",
                        f"{fields.psi[j, i]:.17e}",
                    ]
                )


if __name__ == "__main__":
    main()
