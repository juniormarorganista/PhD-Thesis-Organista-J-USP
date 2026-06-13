"""Generate a parametric ML-ready dataset for MMS Case 12.

Unlike :mod:`03_case12_visualization_and_ml_dataset`, which builds a dataset
for a single ``(Re, Wi, beta, t)`` combination, this example sweeps multiple
combinations and stacks them into a single dataframe, stratified into
train/validation/test splits by parameter combination (``case_id``).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.ml.datasets import (
    build_case12_parametric_dataset,
    save_parametric_dataset_bundle,
)


def _parse_float_list(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(","))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a parametric Case 12 ML dataset (Re, Wi, beta, t sweep)."
    )
    parser.add_argument(
        "--reynolds",
        type=_parse_float_list,
        default=(1.0, 10.0),
        help="Comma-separated Reynolds numbers.",
    )
    parser.add_argument(
        "--wi",
        type=_parse_float_list,
        default=(0.5, 1.0),
        help="Comma-separated Weissenberg numbers.",
    )
    parser.add_argument(
        "--beta",
        type=_parse_float_list,
        default=(0.3, 0.5, 0.7),
        help="Comma-separated beta values.",
    )
    parser.add_argument(
        "--time",
        type=_parse_float_list,
        default=(0.0, 0.5),
        help="Comma-separated time values.",
    )
    parser.add_argument("--nx", type=int, default=17, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=17, help="Number of grid points in y.")
    parser.add_argument("--at", type=float, default=0.0, help="Temporal coefficient.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "case12_parametric",
        help="Output directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    bundle = build_case12_parametric_dataset(
        reynolds_values=args.reynolds,
        wi_values=args.wi,
        beta_values=args.beta,
        time_values=args.time,
        nx=args.nx,
        ny=args.ny,
        at=args.at,
        include_engineered_features=True,
    )
    paths = save_parametric_dataset_bundle(bundle, args.output_dir)

    print("Case 12 parametric ML dataset generated.")
    print(f"Output dir: {args.output_dir}")
    print(f"Combinations: {bundle.metadata['n_combinations']}")
    print(f"Full dataset shape:       {bundle.full.shape}")
    print(f"Train dataset shape:      {bundle.train.shape}")
    print(f"Validation dataset shape: {bundle.validation.shape}")
    print(f"Test dataset shape:       {bundle.test.shape}")
    print(f"Feature columns: {list(bundle.feature_columns)}")
    print(f"Target columns:  {list(bundle.target_columns)}")
    for label, path in paths.items():
        print(f"{label:>10s}: {path}")


if __name__ == "__main__":
    main()
