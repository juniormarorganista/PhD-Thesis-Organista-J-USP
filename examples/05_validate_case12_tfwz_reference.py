"""Validate Python Case 12 `tfwz` against a Fortran reference CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.validation import compare_case12_tfwz_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Python Case 12 tfwz against a Fortran reference CSV."
    )
    parser.add_argument(
        "reference_csv",
        type=Path,
        nargs="?",
        default=Path("validation") / "reference_data" / "case12_tfwz_3x3_default.csv",
    )
    parser.add_argument("--rtol", type=float, default=1.0e-11)
    parser.add_argument("--atol", type=float, default=1.0e-11)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    comparison = compare_case12_tfwz_reference(
        args.reference_csv,
        rtol=args.rtol,
        atol=args.atol,
    )

    print("Case 12 tfwz Fortran-reference validation")
    print(f"reference:      {comparison.reference_path}")
    print(f"grid:           nx={comparison.nx}, ny={comparison.ny}")
    print(f"max_abs_error:  {comparison.max_abs_error:.6e}")
    print(f"mean_abs_error: {comparison.mean_abs_error:.6e}")
    print(f"rtol:           {comparison.rtol:.1e}")
    print(f"atol:           {comparison.atol:.1e}")
    print(f"passed:         {comparison.passed}")

    if not comparison.passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
