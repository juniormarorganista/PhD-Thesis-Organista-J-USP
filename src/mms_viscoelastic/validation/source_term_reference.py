"""Reference-data comparison utilities for manufactured source terms."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import evaluate_case_12_tfwz


@dataclass(frozen=True)
class SourceTermComparison:
    """Summary of a source-term reference comparison."""

    reference_path: Path
    nx: int
    ny: int
    max_abs_error: float
    mean_abs_error: float
    rtol: float
    atol: float
    passed: bool


def load_case12_tfwz_reference(path: str | Path) -> pd.DataFrame:
    """Load a Case 12 `tfwz` reference CSV file."""

    reference_path = Path(path)
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference CSV not found: {reference_path}")

    df = pd.read_csv(reference_path)
    required = {"x", "y", "tfwz"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Reference CSV missing columns: {sorted(missing)}")

    return df.sort_values(["y", "x"]).reset_index(drop=True)


def _infer_grid_shape(reference: pd.DataFrame) -> tuple[int, int]:
    x_values = np.sort(reference["x"].unique())
    y_values = np.sort(reference["y"].unique())
    nx = len(x_values)
    ny = len(y_values)

    if nx * ny != len(reference):
        raise ValueError(
            "Reference CSV does not represent a complete tensor-product grid."
        )

    if not np.isclose(x_values[0], 0.0) or not np.isclose(x_values[-1], 1.0):
        raise ValueError("Only [0, 1] x-domain references are currently supported.")

    if not np.isclose(y_values[0], 0.0) or not np.isclose(y_values[-1], 1.0):
        raise ValueError("Only [0, 1] y-domain references are currently supported.")

    return nx, ny


def compare_case12_tfwz_reference(
    path: str | Path,
    *,
    rtol: float = 1.0e-11,
    atol: float = 1.0e-11,
) -> SourceTermComparison:
    """Compare Python `tfwz` against a Fortran reference CSV."""

    reference_path = Path(path)
    reference = load_case12_tfwz_reference(reference_path)
    nx, ny = _infer_grid_shape(reference)

    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    python_values = evaluate_case_12_tfwz(grid).ravel()
    fortran_values = reference["tfwz"].to_numpy(dtype=float)

    if python_values.shape != fortran_values.shape:
        raise ValueError(
            f"Shape mismatch: Python {python_values.shape}, "
            f"reference {fortran_values.shape}."
        )

    errors = np.abs(python_values - fortran_values)
    max_abs_error = float(np.max(errors))
    mean_abs_error = float(np.mean(errors))
    passed = bool(np.allclose(python_values, fortran_values, rtol=rtol, atol=atol))

    return SourceTermComparison(
        reference_path=reference_path,
        nx=nx,
        ny=ny,
        max_abs_error=max_abs_error,
        mean_abs_error=mean_abs_error,
        rtol=rtol,
        atol=atol,
        passed=passed,
    )
