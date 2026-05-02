"""Generate MMS case-12 figures and ML-ready datasets.

Run from the repository root after installing the package in editable mode:

    python examples/03_case12_visualization_and_ml_dataset.py

Useful options:

    python examples/03_case12_visualization_and_ml_dataset.py --nx 65 --ny 65
    python examples/03_case12_visualization_and_ml_dataset.py --output-dir outputs/case12
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.ml.datasets import build_case12_dataset, save_dataset_bundle
from mms_viscoelastic.visualization import plot_case12_panel, save_case12_component_maps


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate case-12 figures and ML datasets.")
    parser.add_argument("--nx", type=int, default=65, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=65, help="Number of grid points in y.")
    parser.add_argument("--time", type=float, default=0.0, help="Physical time.")
    parser.add_argument("--beta", type=float, default=0.5, help="Solvent-viscosity ratio.")
    parser.add_argument("--at", type=float, default=0.0, help="Temporal exponential coefficient.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "case12",
        help="Directory where CSV files and figures will be written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    parameters = Case12Parameters(beta=args.beta, at=args.at)
    field = evaluate_case_12_exact(grid, t=args.time, parameters=parameters)

    figure_dir = args.output_dir / "figures"
    dataset_dir = args.output_dir / "datasets"

    component_paths = save_case12_component_maps(grid, field, figure_dir)
    panel_path = figure_dir / "case12_panel.png"
    panel = plot_case12_panel(grid, field, output_path=panel_path)
    plt.close(panel)

    table = build_case12_dataset(
        nx=args.nx,
        ny=args.ny,
        time=args.time,
        beta=args.beta,
        at=args.at,
        include_spatial_features=True,
    )
    dataset_paths = save_dataset_bundle(table, dataset_dir, prefix="case12", validation_stride=5)

    print("MMS case 12 visualization + ML dataset")
    print(f"grid: nx={grid.nx}, ny={grid.ny}, points={grid.nx * grid.ny}")
    print(f"output directory: {args.output_dir}")
    print("\nFigures:")
    print(f"  {panel_path}")
    for path in component_paths:
        print(f"  {path}")
    print("\nDatasets:")
    for label, path in dataset_paths.items():
        print(f"  {label:>10s}: {path}")
    print("\nDataset columns:")
    print("  " + ", ".join(table.columns))


if __name__ == "__main__":
    main()
