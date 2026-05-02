"""Example: evaluate MMS case 12 and prepare a tabular ML dataset.

Run from the repository root after installing the package in editable mode:

    python examples/02_case12_manufactured_solution.py

Optional CSV export:

    python examples/02_case12_manufactured_solution.py --export-csv outputs/case12_dataset.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_12_exact
from mms_viscoelastic.ml import field_to_feature_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate MMS case 12 exact fields.")
    parser.add_argument("--nx", type=int, default=33, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=33, help="Number of grid points in y.")
    parser.add_argument("--time", type=float, default=0.0, help="Physical time.")
    parser.add_argument(
        "--export-csv",
        type=Path,
        default=None,
        help="Optional path to export the generated dataset as CSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    field = evaluate_case_12_exact(grid, t=args.time)
    table = field_to_feature_table(grid, field)

    print("MMS case 12 exact field")
    print(f"grid: nx={grid.nx}, ny={grid.ny}, dx={grid.dx:.6e}, dy={grid.dy:.6e}")
    print(table.head())
    print("\ncomponent ranges:")
    for component in field.available_components():
        values = getattr(field, component)
        print(f"  {component:>3s}: min={values.min(): .6e}, max={values.max(): .6e}")

    if args.export_csv is not None:
        args.export_csv.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.export_csv, index=False)
        print(f"\nCSV exported to: {args.export_csv}")


if __name__ == "__main__":
    main()
