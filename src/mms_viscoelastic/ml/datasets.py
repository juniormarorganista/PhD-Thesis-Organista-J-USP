"""Dataset builders for manufactured-solution machine-learning workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.ml.features import TARGET_COLUMNS, field_to_feature_table


DEFAULT_TARGET_COLUMNS = TARGET_COLUMNS


@dataclass(frozen=True)
class DatasetBundle:
    """Container for full/train/validation datasets and metadata."""

    full: pd.DataFrame
    train: pd.DataFrame
    validation: pd.DataFrame
    feature_columns: tuple[str, ...]
    target_columns: tuple[str, ...]
    metadata: dict[str, Any]


def infer_feature_columns(
    df: pd.DataFrame,
    *,
    target_columns: tuple[str, ...] = DEFAULT_TARGET_COLUMNS,
) -> tuple[str, ...]:
    """Return ML input columns by excluding target columns."""

    targets = set(target_columns)
    return tuple(column for column in df.columns if column not in targets)


def split_train_validation(
    df: pd.DataFrame,
    *,
    validation_fraction: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a deterministic train/validation split."""

    if not 0.0 < validation_fraction < 1.0:
        raise ValueError("validation_fraction must be between 0 and 1.")

    validation = df.sample(frac=validation_fraction, random_state=random_state)
    train = df.drop(index=validation.index)

    return (
        train.reset_index(drop=True),
        validation.reset_index(drop=True),
    )


def build_case12_dataset(
    *,
    nx: int = 65,
    ny: int = 65,
    time: float = 0.0,
    beta: float = 0.5,
    at: float = 0.0,
    validation_fraction: float = 0.2,
    random_state: int = 42,
    include_engineered_features: bool = True,
) -> DatasetBundle:
    """Build an ML-ready dataset from the MMS Case 12 exact solution."""

    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    parameters = Case12Parameters(beta=beta, at=at)
    fields = evaluate_case_12_exact(grid, t=time, parameters=parameters)

    full = field_to_feature_table(
        grid,
        fields,
        time=time,
        beta=beta,
        include_engineered_features=include_engineered_features,
    )

    target_columns = tuple(
        column for column in DEFAULT_TARGET_COLUMNS if column in full.columns
    )
    feature_columns = infer_feature_columns(full, target_columns=target_columns)

    train, validation = split_train_validation(
        full,
        validation_fraction=validation_fraction,
        random_state=random_state,
    )

    metadata = {
        "case": "case_12",
        "nx": nx,
        "ny": ny,
        "n_points": len(full),
        "time": time,
        "beta": beta,
        "at": at,
        "validation_fraction": validation_fraction,
        "random_state": random_state,
        "feature_columns": list(feature_columns),
        "target_columns": list(target_columns),
        "shape_convention": "(ny, nx)",
    }

    return DatasetBundle(
        full=full,
        train=train,
        validation=validation,
        feature_columns=feature_columns,
        target_columns=target_columns,
        metadata=metadata,
    )


def save_dataset_bundle(bundle: DatasetBundle, output_dir: str | Path) -> dict[str, Path]:
    """Save dataset CSV files and metadata text file."""

    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)

    full_path = path / "case12_full.csv"
    train_path = path / "case12_train.csv"
    validation_path = path / "case12_validation.csv"
    metadata_path = path / "case12_metadata.txt"

    bundle.full.to_csv(full_path, index=False)
    bundle.train.to_csv(train_path, index=False)
    bundle.validation.to_csv(validation_path, index=False)

    lines = ["Case 12 ML dataset metadata", ""]
    for key, value in bundle.metadata.items():
        lines.append(f"{key}: {value}")
    lines.append("")
    lines.append(f"feature_columns: {list(bundle.feature_columns)}")
    lines.append(f"target_columns: {list(bundle.target_columns)}")

    metadata_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "full": full_path,
        "train": train_path,
        "validation": validation_path,
        "metadata": metadata_path,
    }


__all__ = [
    "DEFAULT_TARGET_COLUMNS",
    "DatasetBundle",
    "build_case12_dataset",
    "infer_feature_columns",
    "save_dataset_bundle",
    "split_train_validation",
]
