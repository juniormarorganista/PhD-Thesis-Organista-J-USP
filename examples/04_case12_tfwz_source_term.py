"""Evaluate the ported MMS Case 12 source terms ``tfwz`` and ``tfTxx``."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import (
    Case12SourceTermParameters,
    evaluate_case_12_tfTxx,
    evaluate_case_12_tfwz,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate Case 12 ported source terms.")
    parser.add_argument("--nx", type=int, default=65)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--time", type=float, default=0.0)
    parser.add_argument("--reynolds", type=float, default=1.0)
    parser.add_argument("--wi", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=0.5)
    parser.add_argument("--epsilon", type=float, default=0.0)
    parser.add_argument("--xi", type=float, default=0.0)
    parser.add_argument("--alpha-g", type=float, default=0.0)
    parser.add_argument("--at", type=float, default=0.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs") / "case12" / "source_terms" / "case12_ported_sources.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    params = Case12SourceTermParameters(
        reynolds=args.reynolds,
        wi=args.wi,
        beta=args.beta,
        epsilon=args.epsilon,
        xi=args.xi,
        alpha_g=args.alpha_g,
        at=args.at,
    )
    tfwz = evaluate_case_12_tfwz(grid, t=args.time, parameters=params)
    tfTxx = evaluate_case_12_tfTxx(grid, t=args.time, parameters=params)

    df = pd.DataFrame(
        {
            "x": grid.X.ravel(),
            "y": grid.Y.ravel(),
            "tfwz": tfwz.ravel(),
            "tfTxx": tfTxx.ravel(),
        }
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)

    print("Case 12 ported source terms generated.")
    print(f"grid: nx={args.nx}, ny={args.ny}")
    for name, array in {"tfwz": tfwz, "tfTxx": tfTxx}.items():
        print(f"{name}:")
        print(f"  shape: {array.shape}")
        print(f"  min:   {np.min(array): .6e}")
        print(f"  max:   {np.max(array): .6e}")
        print(f"  mean:  {np.mean(array): .6e}")
    print(f"csv: {args.output}")


if __name__ == "__main__":
    main()
