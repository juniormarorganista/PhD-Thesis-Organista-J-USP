"""Reference-data comparison utilities for manufactured source terms."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.source_terms import (
    Case12SourceTermParameters,
    evaluate_case_12_tfTxx,
    evaluate_case_12_tfwz,
)


@dataclass(frozen=True)
class SourceTermReferenceComparison:
    """Summary of a numerical comparison against reference data."""

    component: str
    reference_path: Path
    n_points: int
    max_abs_error: float
    max_rel_error: float
    atol: float
    rtol: float
    passed: bool


def load_source_term_reference_csv(path: str | Path) -> pd.DataFrame:
    """Load a source-term reference CSV file."""

    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Reference CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    required = {"x", "y"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Reference CSV is missing required columns: {sorted(missing)}")

    return df


def _grid_from_reference(df: pd.DataFrame) -> Grid2D:
    x_unique = np.sort(df["x"].unique())
    y_unique = np.sort(df["y"].unique())

    nx = len(x_unique)
    ny = len(y_unique)

    if nx * ny != len(df):
        raise ValueError(
            "Reference data must describe a complete structured grid. "
            f"Got nx={nx}, ny={ny}, rows={len(df)}."
        )

    return Grid2D.from_bounds(
        nx=nx,
        ny=ny,
        x0=float(x_unique.min()),
        x1=float(x_unique.max()),
        y0=float(y_unique.min()),
        y1=float(y_unique.max()),
    )


def _relative_error(abs_error: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denominator = np.maximum(1.0, np.abs(reference))
    return abs_error / denominator


def compare_case12_reference_csv(
    path: str | Path,
    *,
    component: str,
    evaluator: Callable[[Grid2D], np.ndarray],
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> SourceTermReferenceComparison:
    """Compare one Case 12 source-term component against a reference CSV."""

    reference_path = Path(path)
    df = load_source_term_reference_csv(reference_path)

    if component not in df.columns:
        raise ValueError(
            f"Reference CSV does not contain component {component!r}. "
            f"Available columns: {list(df.columns)}"
        )

    grid = _grid_from_reference(df)
    predicted = evaluator(grid).ravel()
    reference = df[component].to_numpy(dtype=float)

    if predicted.shape != reference.shape:
        raise ValueError(
            f"Predicted shape {predicted.shape} does not match reference shape {reference.shape}."
        )

    abs_error = np.abs(predicted - reference)
    rel_error = _relative_error(abs_error, reference)

    max_abs = float(abs_error.max(initial=0.0))
    max_rel = float(rel_error.max(initial=0.0))
    passed = bool(np.allclose(predicted, reference, atol=atol, rtol=rtol))

    return SourceTermReferenceComparison(
        component=component,
        reference_path=reference_path,
        n_points=len(df),
        max_abs_error=max_abs,
        max_rel_error=max_rel,
        atol=atol,
        rtol=rtol,
        passed=passed,
    )


def compare_case12_tftxx_reference(
    path: str | Path,
    *,
    parameters: Case12SourceTermParameters | None = None,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> SourceTermReferenceComparison:
    """Compare Python ``tfTxx`` against a Fortran-generated CSV."""

    params = parameters or Case12SourceTermParameters()

    return compare_case12_reference_csv(
        path,
        component="tfTxx",
        evaluator=lambda grid: evaluate_case_12_tfTxx(grid, parameters=params),
        atol=atol,
        rtol=rtol,
    )


def compare_case12_tfwz_reference(
    path: str | Path,
    *,
    parameters: Case12SourceTermParameters | None = None,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> SourceTermReferenceComparison:
    """Compare Python ``tfwz`` against a Fortran-generated CSV."""

    params = parameters or Case12SourceTermParameters()

    return compare_case12_reference_csv(
        path,
        component="tfwz",
        evaluator=lambda grid: evaluate_case_12_tfwz(grid, parameters=params),
        atol=atol,
        rtol=rtol,
    )


__all__ = [
    "SourceTermReferenceComparison",
    "compare_case12_reference_csv",
    "compare_case12_tftxx_reference",
    "compare_case12_tfwz_reference",
    "load_source_term_reference_csv",
]
