"""Generate visualizations and ML-ready datasets for MMS Case 12."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.ml.datasets import build_case12_dataset, save_dataset_bundle
from mms_viscoelastic.visualization.fields import plot_field, plot_field_panel


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Case 12 visualizations and ML datasets."
    )
    parser.add_argument("--nx", type=int, default=65, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=65, help="Number of grid points in y.")
    parser.add_argument("--time", type=float, default=0.0, help="Time value.")
    parser.add_argument("--beta", type=float, default=0.5, help="Model parameter beta.")
    parser.add_argument("--at", type=float, default=0.0, help="Temporal coefficient.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "case12",
        help="Output directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    figures_dir = args.output_dir / "figures"
    datasets_dir = args.output_dir / "datasets"
    figures_dir.mkdir(parents=True, exist_ok=True)
    datasets_dir.mkdir(parents=True, exist_ok=True)

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    parameters = Case12Parameters(beta=args.beta, at=args.at)
    fields = evaluate_case_12_exact(grid, t=args.time, parameters=parameters)

    panel = plot_field_panel(
        grid,
        fields,
        output_path=figures_dir / "case12_panel.png",
        title=f"MMS Case 12 exact fields (nx={args.nx}, ny={args.ny})",
    )
    plt.close(panel)

    for name in ("u", "v", "psi", "wz", "txx", "txy", "tyy"):
        fig, _ = plot_field(
            grid,
            getattr(fields, name),
            output_path=figures_dir / f"case12_{name}.png",
            title=f"Case 12 — {name}",
            label=name,
        )
        plt.close(fig)

    bundle = build_case12_dataset(
        nx=args.nx,
        ny=args.ny,
        time=args.time,
        beta=args.beta,
        at=args.at,
        include_engineered_features=True,
    )
    paths = save_dataset_bundle(bundle, datasets_dir)

    print("Case 12 visualization and ML dataset generated.")
    print(f"Figures:  {figures_dir}")
    print(f"Datasets: {datasets_dir}")
    print(f"Full dataset shape: {bundle.full.shape}")
    print(f"Feature columns: {list(bundle.feature_columns)}")
    print(f"Target columns:  {list(bundle.target_columns)}")
    for label, path in paths.items():
        print(f"{label:>10s}: {path}")


if __name__ == "__main__":
    main()
