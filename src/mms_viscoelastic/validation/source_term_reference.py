"""Reference-data comparison utilities for manufactured source terms.

This module compares Python source-term implementations against CSV files
generated from isolated Fortran reference programs.

Compatibility notes
-------------------
1. The project ``Grid2D.from_bounds`` API currently accepts only ``nx`` and
   ``ny``. Therefore reference CSV coordinates are validated as uniform grids on
   ``[0, 1]`` before building the grid with ``Grid2D.from_bounds(nx, ny)``.

2. Earlier tests and examples expect comparison objects to expose ``nx`` and
   ``ny``. The dataclass keeps these fields explicitly.

3. The original ``tfwz`` validation example also prints ``mean_abs_error``.
   The dataclass therefore exposes both maximum and mean error metrics.
"""

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
    evaluate_case_12_tfTxy,
    evaluate_case_12_tfwz,
)


@dataclass(frozen=True)
class SourceTermReferenceComparison:
    """Summary of a numerical comparison against reference data."""

    component: str
    reference_path: Path
    nx: int
    ny: int
    n_points: int
    max_abs_error: float
    mean_abs_error: float
    max_rel_error: float
    mean_rel_error: float
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


def _assert_uniform_unit_grid(values: np.ndarray, name: str, *, atol: float = 1.0e-12) -> None:
    """Validate that coordinate values describe a uniform grid on [0, 1]."""

    if values.ndim != 1 or len(values) < 1:
        raise ValueError(f"{name} coordinates must be a one-dimensional array.")

    if not np.isclose(values[0], 0.0, atol=atol):
        raise ValueError(
            f"{name} reference coordinates must start at 0.0 for the current "
            f"Grid2D.from_bounds API. Got {values[0]!r}."
        )

    if not np.isclose(values[-1], 1.0, atol=atol):
        raise ValueError(
            f"{name} reference coordinates must end at 1.0 for the current "
            f"Grid2D.from_bounds API. Got {values[-1]!r}."
        )

    if len(values) > 1:
        expected = np.linspace(0.0, 1.0, len(values))
        if not np.allclose(values, expected, atol=atol, rtol=0.0):
            raise ValueError(
                f"{name} reference coordinates must be uniformly spaced on [0, 1]."
            )


def _grid_from_reference(df: pd.DataFrame) -> tuple[Grid2D, int, int]:
    """Build a ``Grid2D`` instance from reference CSV coordinates."""

    x_unique = np.sort(df["x"].unique())
    y_unique = np.sort(df["y"].unique())

    nx = len(x_unique)
    ny = len(y_unique)

    if nx * ny != len(df):
        raise ValueError(
            "Reference data must describe a complete structured grid. "
            f"Got nx={nx}, ny={ny}, rows={len(df)}."
        )

    _assert_uniform_unit_grid(x_unique, "x")
    _assert_uniform_unit_grid(y_unique, "y")

    return Grid2D.from_bounds(nx=nx, ny=ny), nx, ny


def _relative_error(abs_error: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Return a safe relative error with denominator bounded below by one."""

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

    grid, nx, ny = _grid_from_reference(df)
    predicted = evaluator(grid).ravel()
    reference = df[component].to_numpy(dtype=float)

    if predicted.shape != reference.shape:
        raise ValueError(
            f"Predicted shape {predicted.shape} does not match reference shape "
            f"{reference.shape}."
        )

    abs_error = np.abs(predicted - reference)
    rel_error = _relative_error(abs_error, reference)

    max_abs = float(abs_error.max(initial=0.0))
    mean_abs = float(abs_error.mean()) if abs_error.size else 0.0
    max_rel = float(rel_error.max(initial=0.0))
    mean_rel = float(rel_error.mean()) if rel_error.size else 0.0
    passed = bool(np.allclose(predicted, reference, atol=atol, rtol=rtol))

    return SourceTermReferenceComparison(
        component=component,
        reference_path=reference_path,
        nx=nx,
        ny=ny,
        n_points=len(df),
        max_abs_error=max_abs,
        mean_abs_error=mean_abs,
        max_rel_error=max_rel,
        mean_rel_error=mean_rel,
        atol=atol,
        rtol=rtol,
        passed=passed,
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


def compare_case12_tftxy_reference(
    path: str | Path,
    *,
    parameters: Case12SourceTermParameters | None = None,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-10,
) -> SourceTermReferenceComparison:
    """Compare Python ``tfTxy`` against a Fortran-generated CSV."""

    params = parameters or Case12SourceTermParameters()

    return compare_case12_reference_csv(
        path,
        component="tfTxy",
        evaluator=lambda grid: evaluate_case_12_tfTxy(grid, parameters=params),
        atol=atol,
        rtol=rtol,
    )


__all__ = [
    "SourceTermReferenceComparison",
    "compare_case12_reference_csv",
    "compare_case12_tftxx_reference",
    "compare_case12_tftxy_reference",
    "compare_case12_tfwz_reference",
    "load_source_term_reference_csv",
]
