"""Compare Python Case 12 t=0 exact fields against the Fortran reference.

This example validates ``evaluate_case_12_exact`` (the manufactured fields
``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and ``tyy``) against an
isolated Fortran reference program at the initial instant ``t=0``.

Scope note
----------
This is the t=0 initial-instant field comparison described in
``docs/case-12-solver-comparison.md``. A full time-marching solver-level
comparison is left for a future phase.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.validation import CASE12_FIELD_COMPONENTS, compare_case12_field_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate Python Case 12 t=0 exact fields against a Fortran-generated CSV."
    )
    parser.add_argument(
        "reference_csv",
        nargs="?",
        default=Path("validation") / "reference_data" / "case12_fields_3x3_default.csv",
        type=Path,
        help="Reference CSV with columns x,y,u,v,wz,psi,txx,txy,tyy.",
    )
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("Case 12 t=0 exact field reference validation")
    print(f"reference: {args.reference_csv}")
    print()

    all_passed = True
    for component in CASE12_FIELD_COMPONENTS:
        result = compare_case12_field_reference(
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
