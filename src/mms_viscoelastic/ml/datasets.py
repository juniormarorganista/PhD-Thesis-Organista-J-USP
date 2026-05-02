"""Dataset builders for future Machine Learning workflows.

The functions in this module intentionally avoid depending on scikit-learn.
They produce simple, reproducible tabular data that can later feed regression,
operator-learning or physics-informed workflows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.ml.features import field_to_feature_table


def add_spatial_features(table: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic spatial features useful for baseline ML models.

    Added columns:
    - ``x2``, ``y2`` and ``xy``;
    - low-order trigonometric features using ``pi*x`` and ``pi*y``.
    """
    result = table.copy()
    x = result["x"].to_numpy()
    y = result["y"].to_numpy()
    result["x2"] = x**2
    result["y2"] = y**2
    result["xy"] = x * y
    result["sin_pi_x"] = np.sin(np.pi * x)
    result["cos_pi_x"] = np.cos(np.pi * x)
    result["sin_pi_y"] = np.sin(np.pi * y)
    result["cos_pi_y"] = np.cos(np.pi * y)
    return result


def build_case12_dataset(
    *,
    nx: int = 65,
    ny: int = 65,
    time: float = 0.0,
    beta: float = 0.5,
    at: float = 0.0,
    include_spatial_features: bool = True,
) -> pd.DataFrame:
    """Build a tabular dataset from the exact MMS case-12 fields."""
    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    params = Case12Parameters(beta=beta, at=at)
    field = evaluate_case_12_exact(grid, t=time, parameters=params)
    table = field_to_feature_table(grid, field)
    table.insert(2, "time", float(time))
    table.insert(3, "beta", float(beta))
    if include_spatial_features:
        table = add_spatial_features(table)
    return table


def deterministic_spatial_split(
    table: pd.DataFrame,
    *,
    validation_stride: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a table into train/validation sets using grid-index stride.

    This is deterministic and spatially distributed: every ``validation_stride``
    row goes to validation.  It is intentionally simple and reproducible.
    """
    if validation_stride < 2:
        raise ValueError("validation_stride must be at least 2.")
    index = np.arange(len(table))
    validation_mask = index % validation_stride == 0
    validation = table.loc[validation_mask].reset_index(drop=True)
    train = table.loc[~validation_mask].reset_index(drop=True)
    return train, validation


def save_dataset_bundle(
    table: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "case12",
    validation_stride: int = 5,
) -> dict[str, Path]:
    """Save full/train/validation CSV files and a small metadata file."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    train, validation = deterministic_spatial_split(table, validation_stride=validation_stride)

    full_path = output / f"{prefix}_full.csv"
    train_path = output / f"{prefix}_train.csv"
    validation_path = output / f"{prefix}_validation.csv"
    metadata_path = output / f"{prefix}_metadata.txt"

    table.to_csv(full_path, index=False)
    train.to_csv(train_path, index=False)
    validation.to_csv(validation_path, index=False)

    target_columns = ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
    feature_columns = [name for name in table.columns if name not in target_columns]
    metadata = [
        "MMS case 12 dataset bundle",
        f"rows_full={len(table)}",
        f"rows_train={len(train)}",
        f"rows_validation={len(validation)}",
        f"validation_stride={validation_stride}",
        "feature_columns=" + ",".join(feature_columns),
        "target_columns=" + ",".join(target_columns),
    ]
    metadata_path.write_text("\n".join(metadata) + "\n", encoding="utf-8", newline="\n")

    return {
        "full": full_path,
        "train": train_path,
        "validation": validation_path,
        "metadata": metadata_path,
    }
