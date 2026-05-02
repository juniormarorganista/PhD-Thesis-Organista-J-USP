"""Feature-table utilities for machine-learning workflows.

This module contains low-level conversion functions only. It does not import
from :mod:`mms_viscoelastic.ml`, which prevents circular imports.

Default behavior is intentionally conservative: ``field_to_feature_table``
returns only coordinates plus field components unless engineered features are
explicitly requested. This preserves the original test expectation:

``x, y, u, v, psi, wz, txx, txy, tyy``.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd


TARGET_COLUMNS = ("u", "v", "psi", "wz", "txx", "txy", "tyy")


def _as_array(value: Any) -> np.ndarray:
    """Convert *value* to a floating NumPy array."""

    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        raise ValueError("Expected an array-like object, got a scalar.")
    return array


def grid_coordinates(grid: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return ``X, Y`` coordinate arrays using the ``(ny, nx)`` convention."""

    if isinstance(grid, Mapping):
        if "X" in grid and "Y" in grid:
            return _as_array(grid["X"]), _as_array(grid["Y"])
        if "xx" in grid and "yy" in grid:
            return _as_array(grid["xx"]), _as_array(grid["yy"])
        if "x" in grid and "y" in grid:
            x = _as_array(grid["x"])
            y = _as_array(grid["y"])
            if x.ndim == 1 and y.ndim == 1:
                return np.meshgrid(x, y, indexing="xy")
            return x, y

    for x_name, y_name in (("X", "Y"), ("xx", "yy")):
        if hasattr(grid, x_name) and hasattr(grid, y_name):
            return _as_array(getattr(grid, x_name)), _as_array(getattr(grid, y_name))

    if hasattr(grid, "x") and hasattr(grid, "y"):
        x = _as_array(getattr(grid, "x"))
        y = _as_array(getattr(grid, "y"))
        if x.ndim == 1 and y.ndim == 1:
            return np.meshgrid(x, y, indexing="xy")
        return x, y

    if hasattr(grid, "coordinates"):
        coords = grid.coordinates
        if callable(coords):
            coords = coords()
        if isinstance(coords, tuple) and len(coords) == 2:
            return _as_array(coords[0]), _as_array(coords[1])

    raise TypeError(
        "Could not infer grid coordinates. Expected attributes/keys "
        "X/Y, xx/yy, x/y, or a coordinates() method."
    )


def field_mapping(fields: Any) -> dict[str, np.ndarray]:
    """Normalize a field container to a dictionary of 2D arrays."""

    if isinstance(fields, Mapping):
        raw = dict(fields)
    elif is_dataclass(fields):
        raw = asdict(fields)
    else:
        raw = {
            name: getattr(fields, name)
            for name in TARGET_COLUMNS
            if hasattr(fields, name)
        }

    out: dict[str, np.ndarray] = {}
    for name, value in raw.items():
        if name.startswith("_") or value is None:
            continue
        try:
            array = _as_array(value)
        except (TypeError, ValueError):
            continue
        if array.ndim == 2:
            out[name] = array

    if not out:
        raise ValueError("No two-dimensional field arrays were found.")

    return out


def add_engineered_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic coordinate features useful for baseline ML models."""

    out = df.copy()

    if "x" in out.columns:
        out["x2"] = out["x"] ** 2
        out["sin_pi_x"] = np.sin(np.pi * out["x"])
        out["cos_pi_x"] = np.cos(np.pi * out["x"])

    if "y" in out.columns:
        out["y2"] = out["y"] ** 2
        out["sin_pi_y"] = np.sin(np.pi * out["y"])
        out["cos_pi_y"] = np.cos(np.pi * out["y"])

    if "x" in out.columns and "y" in out.columns:
        out["xy"] = out["x"] * out["y"]

    preferred = [
        "x",
        "y",
        "time",
        "beta",
        "x2",
        "y2",
        "xy",
        "sin_pi_x",
        "cos_pi_x",
        "sin_pi_y",
        "cos_pi_y",
    ]
    first = [name for name in preferred if name in out.columns]
    rest = [name for name in out.columns if name not in first]
    return out[first + rest]


def field_to_feature_table(
    grid: Any,
    fields: Any,
    *,
    time: float | None = None,
    beta: float | None = None,
    include_engineered_features: bool = False,
) -> pd.DataFrame:
    """Convert structured fields into a pointwise DataFrame.

    By default, this function returns only:

    ``x, y`` plus the available field components.

    Set ``include_engineered_features=True`` to add polynomial/trigonometric
    coordinate features for machine-learning experiments.
    """

    X, Y = grid_coordinates(grid)
    if X.shape != Y.shape:
        raise ValueError(f"X and Y shapes differ: {X.shape} != {Y.shape}")

    arrays = field_mapping(fields)

    data: dict[str, np.ndarray] = {
        "x": X.ravel(),
        "y": Y.ravel(),
    }

    if time is not None:
        data["time"] = np.full(X.size, float(time))

    if beta is not None:
        data["beta"] = np.full(X.size, float(beta))

    for name in TARGET_COLUMNS:
        if name not in arrays:
            continue
        array = arrays[name]
        if array.shape != X.shape:
            raise ValueError(
                f"Field '{name}' has shape {array.shape}, expected {X.shape}. "
                "Fields must use the NumPy convention (ny, nx)."
            )
        data[name] = array.ravel()

    for name, array in arrays.items():
        if name in TARGET_COLUMNS:
            continue
        if array.shape != X.shape:
            raise ValueError(
                f"Field '{name}' has shape {array.shape}, expected {X.shape}."
            )
        data[name] = array.ravel()

    df = pd.DataFrame(data)

    preferred = ["x", "y", "time", "beta", *TARGET_COLUMNS]
    first = [name for name in preferred if name in df.columns]
    rest = [name for name in df.columns if name not in first]
    df = df[first + rest]

    if include_engineered_features:
        df = add_engineered_features(df)

    return df


__all__ = [
    "TARGET_COLUMNS",
    "add_engineered_features",
    "field_mapping",
    "field_to_feature_table",
    "grid_coordinates",
]
