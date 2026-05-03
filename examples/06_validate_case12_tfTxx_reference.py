"""Validate the Python Case 12 ``tfTxx`` port against Fortran reference data."""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.validation import compare_case12_tftxx_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate Python tfTxx against a Fortran-generated CSV."
    )
    parser.add_argument(
        "reference_csv",
        nargs="?",
        default=Path("validation") / "reference_data" / "case12_tfTxx_3x3_default.csv",
        type=Path,
        help="Reference CSV with columns x,y,tfTxx.",
    )
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = compare_case12_tftxx_reference(
        args.reference_csv,
        atol=args.atol,
        rtol=args.rtol,
    )

    print("Case 12 tfTxx reference validation")
    print(f"reference:     {result.reference_path}")
    print(f"component:     {result.component}")
    print(f"n_points:      {result.n_points}")
    print(f"max_abs_error: {result.max_abs_error:.6e}")
    print(f"max_rel_error: {result.max_rel_error:.6e}")
    print(f"atol:          {result.atol:.6e}")
    print(f"rtol:          {result.rtol:.6e}")
    print(f"passed:        {result.passed}")

    if not result.passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
