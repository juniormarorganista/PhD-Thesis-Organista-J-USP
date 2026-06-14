"""Validate one full Case 12 time-marching step against a Fortran reference CSV.

This example validates ``case12_time_step`` (the explicit Euler step from
Etapa 9, followed by boundary-condition re-imposition, the Poisson re-solve
from Etapa 8 and the recomputation of ``u``/``v`` from the updated ``psi``)
against an isolated Fortran reference program at ``t=dt`` starting from the
``t=0`` exact fields.

Scope note
----------
Both the Python implementation and the Fortran reference use the same
second-order finite-difference operators. The ``wz`` boundary condition is
re-imposed from the exact manufactured solution at ``t+dt`` rather than the
4th-order one-sided stencil used by the original Fortran
``boundarie_vorticity_cavity_flow``. See
``docs/case-12-time-marching-comparison.md`` for the full scope statement.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.validation import (
    CASE12_TIME_STEP_COMPONENTS,
    compare_case12_time_step_reference,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the Python Case 12 time-marching step against a Fortran-generated CSV."
    )
    parser.add_argument(
        "reference_csv",
        nargs="?",
        default=Path("validation") / "reference_data" / "case12_time_step_output_9x9_dt1e-05.csv",
        type=Path,
        help="Reference CSV with columns x,y,wz,txx,txy,tyy,psi,u,v (post-step values).",
    )
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("Case 12 time-marching step Fortran-reference validation")
    print(f"reference: {args.reference_csv}")
    print()

    all_passed = True
    for component in CASE12_TIME_STEP_COMPONENTS:
        result = compare_case12_time_step_reference(
            args.reference_csv,
            component=component,
            atol=args.atol,
            rtol=args.rtol,
        )
        all_passed &= result.passed
        print(
            f"{component:>4s}: max_abs_error={result.max_abs_error:.6e} "
            f"max_rel_error={result.max_rel_error:.6e} passed={result.passed}"
        )

    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
