"""Evaluate MMS Case 01 exact fields and source terms, and validate them
against the Fortran reference CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.io.reference import read_reference_csv
from mms_viscoelastic.manufactured_solutions.case_01 import (
    Case01Parameters,
    evaluate_case_01_exact,
)
from mms_viscoelastic.source_terms.case_01 import (
    Case01SourceTermParameters,
    evaluate_case_01_source_terms,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate Case 01 fields/source terms and compare against a Fortran reference CSV."
    )
    parser.add_argument(
        "reference_csv",
        type=Path,
        nargs="?",
        default=Path("validation") / "reference_data" / "case01_fields_9x9_default.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = read_reference_csv(args.reference_csv)

    nx = len(np.sort(df["x"].unique()))
    ny = len(np.sort(df["y"].unique()))

    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    field = evaluate_case_01_exact(grid, t=0.0, parameters=Case01Parameters())
    terms = evaluate_case_01_source_terms(grid, t=0.0, parameters=Case01SourceTermParameters())

    print("MMS Case 01 - Fortran-reference validation")
    print(f"reference: {args.reference_csv}")
    print(f"grid:      nx={nx}, ny={ny}")
    print()

    columns = {
        "u": field.u,
        "v": field.v,
        "wz": field.wz,
        "psi": field.psi,
        "txx": field.txx,
        "txy": field.txy,
        "tyy": field.tyy,
        "tfwz": terms.tfwz,
        "tfTxx": terms.tfTxx,
        "tfTxy": terms.tfTxy,
        "tfTyy": terms.tfTyy,
    }

    worst = 0.0
    for name, predicted in columns.items():
        reference = df[name].to_numpy(dtype=float).reshape(ny, nx)
        max_abs_error = float(np.max(np.abs(predicted - reference)))
        worst = max(worst, max_abs_error)
        print(f"{name:>6s}: max_abs_error = {max_abs_error:.3e}")

    print()
    print(f"passed: {worst < 1.0e-10}")

    if worst >= 1.0e-10:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
