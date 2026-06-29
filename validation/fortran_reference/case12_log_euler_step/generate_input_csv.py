"""Generate the input CSV consumed by the Case 12 log-Euler-step Fortran reference.

Identical in scope and schema to
``validation/fortran_reference/case12_euler_step/generate_input_csv.py``: it
evaluates the already-validated Python Case 12 exact fields and source
terms at ``t=0`` on a uniform grid and writes::

    x,y,u,v,wz,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy

The log-conformation Fortran reference
(``case12_log_euler_step_reference.f90``) converts ``(txx, txy, tyy)`` to the
log-conformation representation before marching, then converts the
post-step ``Psi`` back to ``(Txx, Txy, Tyy)`` for output, so the same input
schema as the non-log reference is reused.

Usage::

    python validation/fortran_reference/case12_log_euler_step/generate_input_csv.py \\
        validation/reference_data/case12_log_euler_step_input_9x9.csv 9 9
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.source_terms import Case12SourceTermParameters, evaluate_case_12_source_terms

DEFAULT_OUTPUT = "validation/reference_data/case12_log_euler_step_input_9x9.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", default=DEFAULT_OUTPUT)
    parser.add_argument("nx", nargs="?", type=int, default=9)
    parser.add_argument("ny", nargs="?", type=int, default=9)
    parser.add_argument(
        "t",
        nargs="?",
        type=float,
        default=0.3,
        help="Evaluation time. t=0 makes the Case 12 polymeric stress identically zero "
        "everywhere (a degenerate, repeated-eigenvalue conformation matrix), which "
        "exercises neither the eigendecomposition nor the rotation/kernel terms, so "
        "the default is a later time with nonzero stress.",
    )
    args = parser.parse_args()

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    fields = evaluate_case_12_exact(grid, t=args.t, parameters=Case12Parameters())
    source_terms = evaluate_case_12_source_terms(grid, t=args.t, parameters=Case12SourceTermParameters())

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["x", "y", "u", "v", "wz", "txx", "txy", "tyy", "tfwz", "tftxx", "tftxy", "tftyy"]
        )
        for j in range(grid.ny):
            for i in range(grid.nx):
                writer.writerow(
                    [
                        f"{grid.X[j, i]:.17e}",
                        f"{grid.Y[j, i]:.17e}",
                        f"{fields.u[j, i]:.17e}",
                        f"{fields.v[j, i]:.17e}",
                        f"{fields.wz[j, i]:.17e}",
                        f"{fields.txx[j, i]:.17e}",
                        f"{fields.txy[j, i]:.17e}",
                        f"{fields.tyy[j, i]:.17e}",
                        f"{source_terms.tfwz[j, i]:.17e}",
                        f"{source_terms.tfTxx[j, i]:.17e}",
                        f"{source_terms.tfTxy[j, i]:.17e}",
                        f"{source_terms.tfTyy[j, i]:.17e}",
                    ]
                )


if __name__ == "__main__":
    main()
