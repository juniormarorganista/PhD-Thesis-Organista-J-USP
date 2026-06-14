"""Validate multi-step Case 12 time marching against a Fortran reference CSV.

This example validates ``case12_multistep_march`` (Etapa 11), which chains
``case12_time_step`` (Etapa 10) over ``n_steps`` increments of size ``dt``
starting from the ``t=0`` exact fields, against an isolated Fortran reference
program after the same number of steps.

Scope note
----------
Both the Python implementation and the Fortran reference use the same
second-order finite-difference operators and ``wz``-boundary simplification as
Etapa 10. With the default manufactured-solution parameter ``at=0``, the exact
boundary fields and the manufactured source terms are independent of ``t``, so
the Fortran reference reuses the ``t=0`` input CSV at every step. See
``docs/case-12-multistep-marching-comparison.md`` for the full scope
statement.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.validation import (
    CASE12_TIME_STEP_COMPONENTS,
    compare_case12_multistep_reference,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the Python Case 12 multi-step time marching against a Fortran-generated CSV."
    )
    parser.add_argument(
        "reference_csv",
        nargs="?",
        default=Path("validation") / "reference_data" / "case12_multistep_output_9x9_dt1e-05_n10.csv",
        type=Path,
        help="Reference CSV with columns x,y,wz,txx,txy,tyy,psi,u,v (values after n_steps steps).",
    )
    parser.add_argument("--n-steps", type=int, default=10)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("Case 12 multi-step time-marching Fortran-reference validation")
    print(f"reference: {args.reference_csv}")
    print(f"n_steps:   {args.n_steps}")
    print()

    all_passed = True
    for component in CASE12_TIME_STEP_COMPONENTS:
        result = compare_case12_multistep_reference(
            args.reference_csv,
            n_steps=args.n_steps,
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
